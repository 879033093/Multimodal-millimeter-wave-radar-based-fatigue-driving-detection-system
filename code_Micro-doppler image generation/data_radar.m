function [retVal] = data_radar(adcdata,n_samples,n_chirps,xx)

data_radar_1 = reshape(adcdata(1,:),n_samples,[]);   %RX1
data_radar_2 = reshape(adcdata(2,:),n_samples,[]);   %RX2
data_radar_3 = reshape(adcdata(3,:),n_samples,[]);   %RX3
data_radar_4 = reshape(adcdata(4,:),n_samples,[]);   %RX4
data_radar=zeros(n_samples,n_chirps*n_samples,4);            
data_radar(:,:,1)=data_radar_1;   
data_radar(:,:,2)=data_radar_2;
data_radar(:,:,3)=data_radar_3;
data_radar(:,:,4)=data_radar_4;

data_radar=data_radar(:,(1+128*(xx-1)):(128*(xx-1)+n_chirps),:);
retVal = data_radar;