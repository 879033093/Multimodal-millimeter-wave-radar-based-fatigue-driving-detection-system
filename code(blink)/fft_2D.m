function [retVal] = fft_2D(data_rxn, n_samples, n_chirps, No_frame)
    data = zeros(n_samples, n_chirps);
    range_win = hamming(n_samples);
    doppler_win = hamming(n_chirps);

    data_temp = data_rxn(:,:,No_frame)' - mean((data_rxn(:,:,No_frame)'));
    data_temp = data_temp';

    for j = 1:n_chirps
        temp = data_temp(:,j) .* range_win;
        temp_fft = fft(temp, n_samples);
        data(:,j) = temp_fft;
    end

    ns = size(data,2);
    [b,a] = butter(4, 0.0075, 'high');
    for k = 1:size(data,1)
        data(k,1:ns) = filter(b,a,data(k,1:ns));
    end

    for k = 1:n_samples
        temp = data(k,:) .* (doppler_win)';
        temp_fft = fftshift(fft(temp, n_chirps));
        data(k,:) = temp_fft;
    end

    retVal = data;
end
