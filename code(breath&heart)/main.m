close all;
clc;
clear all;

f0=77e8;
numADCSamples = 200;
Doppler_Number = 200;
freqSlopeRate = 70.006e12;
adcSampleRate = 4e6;

Ts = 50e-3;
Tc = numADCSamples/adcSampleRate;
B = Tc*freqSlopeRate;
deltaR = 3e8/(2*B);

n_chirps=2;
n_samples=200;
N=200;
M=2;
fs=4e6;
c=3.0e8;
f0=77e9;
lambda=c/f0;
d=lambda/2;
Tc=57e-6;
Tf=50e-3;
B=3990.34e6;
S=B/(Tc-(7e-6));
Bv=(N/fs)*S;
rangeRes=c/2/Bv;
Ta=0:Tf:100*Tf;
Rr=0:rangeRes:(N-1)*rangeRes;
Vr=(-M/2:M/2-1)*lambda/Tc/M/2;

numADCSamples = 200;
numADCBits = 16;
numRX = 8;
numTx =2;
numLanes = 2;
isReal = 0;
chirpLoop = 2;

for tianxian=1:8     
    addpath('      ');% path of row data
    Filename = 'huxin5.bin';
    fid = fopen(Filename,'r');
    adcDataRow = fread(fid, 'int16');
    if numADCBits ~= 16
        l_max = 2^(numADCBits-1)-1;
        adcDataRow(adcDataRow > l_max) = adcDataRow(adcDataRow > l_max) - 2^numADCBits;
    end
    fclose(fid);

    fileSize = size(adcDataRow, 1);
    PRTnum = fix(fileSize/(numADCSamples*numRX));
    fileSize = PRTnum * numADCSamples*numRX;
    adcData = adcDataRow(1:fileSize);
    if isReal
        numChirps = fileSize/numADCSamples/numRX;
        LVDS = zeros(1, fileSize);
        LVDS = reshape(adcData, numADCSamples*numRX, numChirps);
        LVDS = LVDS.';
    else
        numChirps = fileSize/2/numADCSamples/numRX;
        LVDS = zeros(1, fileSize/2);
        counter = 1;
        for i=1:4:fileSize-1
            LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2);
            LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3); counter = counter + 2;
        end
        LVDS = reshape(LVDS, numADCSamples*numRX, numChirps);
        LVDS = LVDS.';
    end

    numframe=1024;
    adcData = zeros(numRX,numChirps*numADCSamples);
    for row = 1:numRX
        for i = 1:numChirps
            adcData(row, (i-1)*numADCSamples+1:i*numADCSamples) = LVDS(i, (row-1)*numADCSamples+1:row*numADCSamples);
        end
    end

    adcdata=adcData;

    n_frame=floor(length(adcdata)/n_chirps/n_samples);
    retVal= reshape(adcData(tianxian, :), numADCSamples, numChirps);

    process_adc=zeros(numADCSamples,numChirps/2);
    for nchirp = 1:2:numChirps
        process_adc(:, (nchirp-1)/2+1) = retVal(:,nchirp);
    end
    channel_data(:,:,tianxian)=process_adc';
end

rangeFFTNum=256;
useFramesNum=1024;
loop_cnt = floor(length(channel_data(:,1,1)) / useFramesNum);

k = 1;

use_channel_data = channel_data((k - 1) * useFramesNum + 1 : k * useFramesNum,:,:);

rangeProfile = MTI_PulseCompression(use_channel_data,0,rangeFFTNum);
sum_rangProfile = sum(abs(rangeProfile(:,:,1)),1);
[~,targetIndex] = max(sum_rangProfile);

searchAngleRange = 60;

[~,azimuSpectrogram,Rxv] = IWR1642ODS_DOA(rangeProfile,2,useFramesNum,searchAngleRange);

maxAzimu = max(azimuSpectrogram,[],2);

[values,peaks_index] = findpeaks(maxAzimu,'minpeakheight',500000);

figure(1);  
plot(1:length(maxAzimu),10*log10(maxAzimu));    
hold on; grid on
plot(peaks_index,10*log10(maxAzimu(peaks_index)),'bd');
xlabel('ang/°');ylabel('G（dB）');title('MVDR');
hold off

figure(2);
imagesc(-searchAngleRange:searchAngleRange,(1 :rangeFFTNum)*deltaR,abs(azimuSpectrogram.'));
ylabel(' d(m) ');xlabel(' ang(°)');title('d——ang');
axis xy

[values_sort,index] =sort(values,'descend');
peaks_index_sort =peaks_index(index);
peaks_index_max =peaks_index_sort(1:2);

distance_indices = zeros(1, length(peaks_index_max));
distance_values = zeros(1, length(peaks_index_max));

for dd = 1:length(peaks_index_max)
    [distance_values(dd), distance_indices(dd)] = max(azimuSpectrogram(peaks_index_max(dd), :));
end

target_rangeProfile  = zeros(useFramesNum,length(peaks_index_max));

for m = 1:length(peaks_index_max)

    detAngle = -searchAngleRange + peaks_index_max(m) * (searchAngleRange * 2 / length(azimuSpectrogram(:,1)));
    dashabi111(:,m)=   detAngle;

    fai = 2 * pi * sin(detAngle / 180 * pi) * d / lambda;
    dashabi222(:,m)= fai;
    aTheta = [1,exp(-1j*1*fai),exp(-1j*2*fai),exp(-1j*3*fai),exp(-1j*4*fai),exp(-1j*5*fai),exp(-1j*6*fai),exp(-1j*7*fai)].';

    Wopt = (Rxv  * aTheta) / (aTheta' * Rxv  * aTheta);   

    xt = squeeze(rangeProfile(:,distance_indices(:,m),:));
    target_rangeProfile(:,m) = xt * Wopt; 
end

[breathRate,heartRate] = get_heartBreath_rate(target_rangeProfile,1/Ts);
