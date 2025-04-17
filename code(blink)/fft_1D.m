function [retVal] = fft_1D(data_rxn, n_samples, n_chirps, No_frame)
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

    retVal = data;
end
