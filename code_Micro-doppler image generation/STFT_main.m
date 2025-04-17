
data_rx4_h=reshape(rx4, n_samples, n_chirps* n_frame); 
ReIm_Data=data_rx4_h' ;
jingtai = zeros(M*60,N);            
    avg = sum(ReIm_Data(:,:),2)/200;
    for chirp=1:200
       jingtai(:,chirp) = ReIm_Data(:,chirp)-avg;
    end
jingtai=jingtai';


Data_time=jingtai;     
win = ones(N,7680);
tmp = fftshift(fft(Data_time.*win),1);  
                                        
                                        
                                        
Data_range=zeros(size(tmp,1)/2,size(tmp,2)); 
Data_range(1:N/2,:) = tmp(N/2+1:N,:);
                                     

ns = oddnumber(size(Data_range,2))-1;        
Data_range_MTI = zeros(size(Data_range,1),ns);


PRF=20;
f_low = 0.075;
f_high = 2;
[b, a] = butter(4, [f_low, f_high] / (PRF / 2), 'bandpass');

% 计算频率响应
[h, f1] = freqz(b, a, ns);

%b 和 a 是 IIR 滤波器的分子和分母系数。它们定义了滤波器的传递函数。
%h 是滤波器的频率响应，通常是一个复数向量，表示滤波器在不同频率处的增益和相位信息。
%f1 是频率响应的频率轴。
%ns 是频率响应的计算点数。
for k=1:size(Data_range,1)    %k=100
  Data_range_MTI(k,1:ns) = filter(b,a,Data_range(k,1:ns));   %对列进行滤波（每列的50000个数据）
end

freq =(0:ns-1)*fs/(2*ns); 
%range_axis=(freq*3e8*Tsweep)/(2*Bw);





%% stft阶段
                  
bin_indl = 4;       
bin_indu = 10;


PRF=round(1/Tc);         
TimeWindowLength = 200;   
OverlapFactor = 0.95;    
OverlapLength = round(TimeWindowLength*OverlapFactor); 
Pad_Factor = 4;                   %
FFTPoints = Pad_Factor*TimeWindowLength;
DopplerBin=PRF/(FFTPoints);%计算每个Doppler速度单元的宽度
DopplerAxis=-PRF/2:DopplerBin:PRF/2-DopplerBin;%生成Doppler频率轴。
WholeDuration=size(Data_range_MTI,2)/PRF;%计算整个信号的持续时间。
NumSegments=floor((size(Data_range_MTI,2)-TimeWindowLength)/floor(TimeWindowLength*(1-OverlapFactor)));%计算可用于多普勒谱的时间窗段数，用于后续循环操作
%计算可用于多普勒谱的时间窗段数，考虑了重叠因子。
Data_spec_MTI2=0;
Data_spec=0;

for RBin=bin_indl:1:bin_indu  

    %Data_range_MTI（63*1000）是距离计算中的结果  spectrogram是stft函数输入参数（窗长，重叠长，fft点数）
    Data_MTI_temp = fftshift(spectrogram(Data_range_MTI(RBin,:),TimeWindowLength,OverlapLength,FFTPoints),1);
   %对Data_range_MTI(RBin,:)的列短时傅里叶，从大小1*10000变成800*981
    Data_spec_MTI2=Data_spec_MTI2+abs(Data_MTI_temp);  %重复3-10个距离单元

    %下部分是没经过MIT的
    %Data_range_MTI（63*1000）是距离计算中的结果  spectrogram是stft函数输入参数（窗长，重叠长，fft点数）
    Data_temp = fftshift(spectrogram(Data_range(RBin,:),TimeWindowLength,OverlapLength,FFTPoints),1);
   %对Data_range_MTI(RBin,:)的列短时傅里叶，从大小1*10000变成800*981
    Data_spec=Data_spec+abs(Data_temp);  %重复10-20个距离单元


end


figure
imagesc(Ta,Vr,20*log10(abs(Data_spec_MTI2))); colormap('jet'); axis xy
ylim([-6 6]);
clim = get(gca,'CLim');
set(gca, 'CLim', clim(2)+[-40,0]);
xlabel('Time[s]', 'FontSize',16);
ylabel('Velocity [m/s]','FontSize',16)
set(gca, 'FontSize',16)



