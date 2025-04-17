function [retVal] = readdata(fileName)



numADCBits = 16; 
numLanes = 4; 
isReal = 0; 

fid=fopen(fileName,'rb'); 
adcData = fread(fid, 'int16');

if numADCBits ~= 16
    l_max = 2^(numADCBits-1)-1;
    adcData(adcData > l_max) = adcData(adcData > l_max) - 2^numADCBits;
end

fclose(fid);

if isReal==1

   adcData = reshape(adcData, numLanes, []);

else

adcData = reshape(adcData, numLanes*2, []);
adcData = adcData([1,2,3,4],:) + sqrt(-1)*adcData([5,6,7,8],:);
end

retVal = adcData;

end
