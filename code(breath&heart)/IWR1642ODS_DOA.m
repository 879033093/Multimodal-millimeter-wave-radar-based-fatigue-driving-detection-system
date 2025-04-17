function [pitchSpectrogram, azimuSpectrogram, Rxv] = IWR1642ODS_DOA(rangeProfile, doaFlag, useChirpNum, searchAngleRange, wopt_flag, wopt)

FFT_DOA_MODE   = 0;
CBF_DOA_MODE   = 1;
MVDR_DOA_MODE  = 2;

rangNum = size(rangeProfile, 2);
lambda = 5e-3;
d = lambda / 2;
if (doaFlag == CBF_DOA_MODE) || (doaFlag == MVDR_DOA_MODE)
    SEARCH_ANGLE_RANGE = searchAngleRange / 180 * pi;
    SEARCH_ANGLE_SPACE = 1 / 180 * pi;
    fai = 2 * pi * sin(-SEARCH_ANGLE_RANGE:SEARCH_ANGLE_SPACE:SEARCH_ANGLE_RANGE) * d / lambda;
    azimuSpectrogram = zeros(length(fai), rangNum);
    pitchSpectrogram = zeros(length(fai), rangNum);
    for r = 1:rangNum
        xt = squeeze(rangeProfile(1:useChirpNum, r, :))';
        Rx = xt * xt';
        for an = 1:length(fai)
            aTheta = [1, exp(-j*1*fai(an)), exp(-j*2*fai(an)), exp(-j*3*fai(an)), exp(-j*4*fai(an)), exp(-j*5*fai(an)), exp(-j*6*fai(an)), exp(-j*7*fai(an))];
            if (doaFlag == CBF_DOA_MODE)
                azimuSpectrogram(an, r) = abs(aTheta * Rx * aTheta'); 
                Rxv = pinv(Rx);  
            else
                Rxv = pinv(Rx);  
                azimuSpectrogram(an, r) = 1 / abs(aTheta * Rxv * aTheta');
            end
        end
    end
end
AzimuSpectrogram = azimuSpectrogram / max(azimuSpectrogram);
end
