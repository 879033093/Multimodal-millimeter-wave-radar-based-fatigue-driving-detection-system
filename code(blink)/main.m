 clc; clear;close all;
                                       
%% Parameter Settings
n_chirps=255;                   
n_samples=128;                  
N=128;                         
M=255;                          
fs=2e6;                       
c=3.0e8;                       
f0=77e9;                        
lambda=c/f0;                    
d=lambda/2;                     
Tc=114e-6;                       
Tf=100e-3;                       
B=3990.34e6;                    
S=B/(Tc-(14e-6));
Bv=(N/fs)*S;
rangeRes=c/2/Bv;                
Ts=0:Tc:127*Tc;                 
Ta=0:Tf:30*Tf;                 
Rr=0:rangeRes:(N-1)*rangeRes;   
Vr=(-M/2:M/2-1)*lambda/Tc/M/2;  

%% Read File
fname='wang2(2).bin'; %Read File Names
adcdata =readdata(fname);    

rx1=abs(adcdata(1,:));
rx2=abs(adcdata(2,:));
rx3=abs(adcdata(3,:));
rx4=abs(adcdata(4,:));
RX=rx1+rx2+rx3+rx4;
n_frame=floor(length(adcdata)/n_chirps/n_samples); 
                                
data_rx1= reshape(adcdata(1,:),n_samples,n_chirps,n_frame); 
data_rx2= reshape(adcdata(2,:),n_samples,n_chirps,n_frame); 
data_rx3= reshape(adcdata(3,:),n_samples,n_chirps,n_frame);
data_rx4= reshape(adcdata(4,:),n_samples,n_chirps,n_frame); 

data_RX= reshape(RX,n_samples,n_chirps,n_frame); 
data_RX_h=reshape(data_RX, 128, 255* 30);  
 
range_rx1_data=data_rx2(:,:,20) ;

fft1d= zeros(M,N);

    for chirp_fft=1:M
        fft1d(chirp_fft,:) = db(abs(fft((range_rx1_data(:,chirp_fft)))));
    end

%Micro-Doppler 
profile=zeros(n_frame,n_chirps);    
for n=1:1:30 
   
    rx1_2dfft=fft_2D(data_rx1,n_samples,n_chirps,n); 

    temp=10*log(abs(rx1_2dfft(:,:,1)));        

    rx2_2dfft=fft_2D(data_rx2,n_samples,n_chirps,n);
    rx3_2dfft=fft_2D(data_rx3,n_samples,n_chirps,n);
    rx4_2dfft=fft_2D(data_rx4,n_samples,n_chirps,n);

    rx_RDM(:,:,1)=rx1_2dfft;
    rx_RDM(:,:,2)=rx2_2dfft;
    rx_RDM(:,:,3)=rx3_2dfft;
    rx_RDM(:,:,4)=rx4_2dfft;

    rx_2dfft=abs(rx1_2dfft)+abs(rx2_2dfft)+abs(rx3_2dfft)+abs(rx4_2dfft);
    temp=10*log(abs(rx_2dfft));      
    %2D_CFAR
    rx_2dcfar=CFAR_2D(rx_2dfft);      
    [rx_2dcfar_plots,add,num_target]=PlotsCentroid(rx_2dcfar);
            
   %Micro-Doppler
   for j=2:5    

      profile(n,:)=profile(n,:)+ rx_2dcfar(j,:);  
                                                
   end
                                               
hold off

end
% Micro-Doppler Map
profile=profile';
figure(5)
colormap(jet);
imagesc(Ta,Vr,(profile));    
title('Micro-Doppler Map','FontSize', 14);xlabel('Time s','FontSize', 14);ylabel('Velocity  m/s','FontSize', 14);

%% 


[X,Y] = meshgrid((0:N-1)*rangeRes, ...        
   (1:M*30));                          

fft1d_SUM= zeros(M*30,N);

    for chirp_fft=1:M*30
        fft1d_SUM(chirp_fft,:) = fft((data_RX_h(:,chirp_fft)));%对抽取出来的1024个chirp分别fft
    end

ns = size(fft1d_SUM,2);
[b,a] = butter(4, 0.0075, 'high');
[h, f1] = freqz(b, a, ns);

for k=1:size(fft1d_SUM,1)
  fft1d_SUM(k,1:ns) = filter(b,a,fft1d_SUM(k,1:ns)); 
end
fft1d_SUM(1,:)=0;
fft1d_SUM=fft1d_SUM';
fft1d_SUM_to_look_look=sum(fft1d_SUM(2:7,:));
figure
plot(abs(fft1d_SUM_to_look_look));
