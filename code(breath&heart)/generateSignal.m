function rawData = generateSignal(Parameter)

    c = Parameter.c;                
    stratFreq = Parameter.stratFreq; 

    Tr = Parameter.Tr;            
    samples = Parameter.Samples;  
    fs = Parameter.Fs;        

    rangeBin = Parameter.rangeBin;     
    chirps = Parameter.Chirps;        
    dopplerBin = Parameter.dopplerBin; 

    slope = Parameter.Slope;          
    bandwidth = Parameter.Bandwidth;   
    centerFreq = Parameter.centerFreq; 
    lambda = Parameter.lambda;
    txAntenna = Parameter.txAntenna; 
    txNum = length(txAntenna);       
    rxAntenna = Parameter.rxAntenna;
    rxNum = length(rxAntenna);       
    dz = Parameter.dz;          
    dx = Parameter.dx;          
    
    target = Parameter.target;  
    targetNum = size(target,1); 
    rawData = zeros(txNum*rxNum,rangeBin,dopplerBin);

    t = 0:1/fs:Tr-(1/fs); 
    for chirpId = 1:chirps
       for txId = 1:txNum 
            St = exp((1i*2*pi)*(centerFreq*(t+((txNum-1)*chirps + chirpId)*Tr)+slope/2*t.^2)); 
            for rxId = 1:rxNum
                Sif = zeros(1,rangeBin);
                for targetId = 1:targetNum

                    
                    if targetId==1
                        targetRange = target(targetId,1)-Parameter.frame; 
                        targetSpeed = target(targetId,2); 
                        targetAngle = target(targetId,3);
                    elseif targetId==2
                        targetRange = target(targetId,1)+0.5*Parameter.frame; 
                        targetSpeed = target(targetId,2); 
                        targetAngle = target(targetId,3);
                    elseif targetId==3
                        targetRange = target(targetId,1)+Parameter.frame; 
                        targetSpeed = target(targetId,2); 
                        targetAngle = target(targetId,3);
                    end

                    tau = 2 * (targetRange + targetSpeed * (txId - 1) * Tr) / c;
                    fd = 2 * targetSpeed / lambda;
                    wx = ((txId-1) * rxNum + rxId) / lambda * dx * sind(targetAngle);
                    Sr = 10*exp((1i*2*pi)*((centerFreq-fd)*(t-tau+((txNum-1)*chirps + chirpId) * Tr)+slope/2*(t-tau).^2 + wx));  %回波信号
                    Sif = Sif + St .* conj(Sr);
                    
                    Sif = awgn(Sif,20);
                end
                rawData((txId-1) * rxNum + rxId,:,chirpId) = Sif;
            end
        end
    end
end
