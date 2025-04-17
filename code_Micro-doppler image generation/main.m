clear all
clc

    n_chirps = 255;                   
    n_samples = 200;                  
    n_frame = 60;
    N = 200;                         
    M = 255;                          
    fs = 4e6;                       
    c = 3.0e8;                        
    f0 = 77e9;                        
    lambda = c / f0;                  
    d = lambda / 2;                   
    Tc = 64e-6;                       
    Tf = 50e-3;                       
    B = 3990.3e6;                     
    S = B / (Tc - (7e-6));
    Bv = (N / fs) * S;
    rangeRes = c / 2 / Bv;           
    Ts = 0:Tc:255 * Tc;               
    Ta = 0:Tf:59 * Tf;               
    Rr = 0:rangeRes:(N - 1) * rangeRes; 
    Vr = (-M / 2:M / 2 - 1) * lambda / Tc / M / 2; 
    numADCBits = 16;     
    numRX = 4;           
    numLanes = 2;        
    isReal = 0;          


%% read Bin
Filename = '1.bin';  

 filePath = fullfile('D:\1.研究生\个人数据测量\第二篇\头部动作原始数据集\头部动作bin\车内\点头\何', Filename);
    fileID = fopen(filePath, 'rb'); 
    
    if fileID == -1
        fprintf('???: %s\n', filePath);
       
    end
    
    data = fread(fileID, 'int16'); %

  
    adcDataRow = data;
    if numADCBits ~= 16 && numADCBits <= 16 %   
        l_max = 2^(numADCBits - 1) - 1;  
        offset = 2^(numADCBits - 1);  
        adcDataRow(adcDataRow > l_max) = adcDataRow(adcDataRow > l_max) - offset;   
    end

    fileSize = size(adcDataRow, 1);
    PRTnum = fix(fileSize / (n_samples * numRX));
    fileSize = PRTnum * n_samples * numRX;
    adcData = adcDataRow(1:fileSize);
    
    if isReal
        numChirps = fileSize / n_samples / numRX;
        LVDS = reshape(adcData, n_samples * numRX, numChirps).';
    else
        numChirps = fileSize / 2 / n_samples / numRX; 
        LVDS = zeros(1, fileSize / 2);
        counter = 1;
        for j = 1:4:fileSize - 1
            LVDS(counter) = adcData(j) + 1i * adcData(j + 2);
            LVDS(counter + 1) = adcData(j + 1) + 1i * adcData(j + 3);
            counter = counter + 2;
        end
        LVDS = reshape(LVDS, n_samples * numRX, numChirps).';
    end

   
    adcData = zeros(numRX, numChirps * n_samples);
    for row = 1:numRX
        for j = 1:numChirps
            adcData(row, (j - 1) * n_samples + 1:j * n_samples) = LVDS(j, (row - 1) * n_samples + 1:row * n_samples);
        end
    end

    data_rx1 = reshape(adcData(1, :), n_samples, n_chirps, n_frame); 
    data_rx4 = reshape(adcData(4, :), n_samples, n_chirps, n_frame);

    %% 
    data_rx4_h = reshape(data_rx4, n_samples, n_chirps * n_frame); 
    ReIm_Data = data_rx4_h';
    avg = mean(ReIm_Data, 2);
    jingtai = ReIm_Data - avg;
    jingtai = jingtai';

    %% 
    Data_time = jingtai; 
    win = ones(N, 15300);
    tmp = fftshift(fft(Data_time .* win), 1);

    Data_range = tmp(N/2 + 1:N, :);
    ns = size(Data_range, 2) - 1;
    Data_range_MTI = zeros(size(Data_range, 1), ns);

    PRF = 20;
    f_low = 0.15;
    f_high = 2;
    [b, a] = butter(4, [f_low, f_high] / (PRF / 2), 'bandpass');

    for k = 1:size(Data_range, 1)
        Data_range_MTI(k, 1:ns) = filter(b, a, Data_range(k, 1:ns));
    end

    %% stft阶段
    bin_indl = 4;
    bin_indu = 9;
    TimeWindowLength = 200;
    OverlapFactor = 0.95;
    OverlapLength = round(TimeWindowLength * OverlapFactor);
    Pad_Factor = 4;
    FFTPoints = Pad_Factor * TimeWindowLength;
    DopplerBin = PRF / FFTPoints;
    DopplerAxis = -PRF / 2:DopplerBin:PRF / 2 - DopplerBin;
    WholeDuration = size(Data_range_MTI, 2) / PRF;
    NumSegments = floor((size(Data_range_MTI, 2) - TimeWindowLength) / floor(TimeWindowLength * (1 - OverlapFactor)));
    Data_spec_MTI2 = 0;
    Data_spec = 0;

    for RBin = bin_indl:bin_indu
        Data_MTI_temp = fftshift(spectrogram(Data_range_MTI(RBin, :), TimeWindowLength, OverlapLength, FFTPoints), 1);
        Data_spec_MTI2 = Data_spec_MTI2 + abs(Data_MTI_temp);
    end

    
    figure;
    imagesc(Ta, Vr, 20 * log10(abs(Data_spec_MTI2)));
    colormap('jet');
    axis xy;
    ylim([-6 6]);
    clim = get(gca, 'CLim');
    set(gca, 'CLim', clim(2) + [-40, 0]);
    set(gca, 'XTick', [], 'YTick', []);
    set(gca, 'XColor', 'none', 'YColor', 'none');
