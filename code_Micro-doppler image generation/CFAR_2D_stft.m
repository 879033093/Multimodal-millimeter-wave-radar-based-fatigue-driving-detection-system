function [retVal] = CFAR_2D_stft(RDM)

    Tr = 8; 
    Td = 4;
    Gr = 4;
    Gd = 2;

    [r, c] = size(RDM);
    CFAR = zeros(r, c);

    r_margin = 2 * (Tr + Gr);
    d_margin = 2 * (Td + Gd);
    gridSize = (2 * Tr + 2 * Gr + 1) * (2 * Td + 2 * Gd + 1);
    numGcells = (2 * Gr + 1) * (2 * Gd + 1);
    numTcells = gridSize - numGcells;

    pfa = 1e-10;
    alpha = numTcells * (pfa^(-1 / numTcells) - 1);
    snr_offset = 16 * log10(alpha);

    RDM_temp1 = zeros(r, c);
    for i = 1:r
        for j = 1:c
            RDM_temp1(i, j) = 10 * log10(abs(RDM(i, j)));
        end
    end

    a = mean(RDM_temp1, 'all');
    RDM_temp = padarray(RDM_temp1, [Tr + Gr, Td + Gd], a, 'both');
    noise_level = zeros(r, c);

    for i = 1:r
        for j = 1:c
            sig_pow = (RDM_temp(i:i + r_margin, j:j + d_margin));
            G_pow = (RDM_temp(i + Tr:i + Tr + Gr * 2, j + Td:j + Td + Gd * 2));
            sig_sum = sum(sum(sig_pow)) - sum(sum(G_pow));
            noise_level(i, j) = (sig_sum / numTcells);
            sig_threshold = noise_level(i, j) + snr_offset;

            if (RDM_temp(i + r_margin / 2, j + d_margin / 2) > sig_threshold)
                CFAR(i, j) = RDM_temp1(i, j);
            end  
        end
    end

    retVal = CFAR;
end
