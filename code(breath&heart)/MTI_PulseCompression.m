function Cancel_PulseCompre_Data = MTI_PulseCompression(framesData,cancelFlag,rangeFFTNum,pluseCancelSpace)
% 内容：对一帧雷达回波做对消和脉冲压缩处理

% 输入：   
% framesData：带对消和脉冲压缩的雷达回波数据
% cancelFlag：对消标志位：cancelFlag=0 均值对消；cancelFlag=1 两脉冲对消；cancelFlag=2 三脉冲对消
% rangeFFTNum：脉压FFT点数
% pluseCancelSpace：两脉冲对消间隔，只有用到两脉冲对消时才需要传递此参数

% 输出：
% Cancel_PulseCompre_Data：做了脉压之后的复数数据

%% 对数据进行对消和脉冲压缩
meanCancel_mode = 0;
twopCancel_mode = 1;
thrpCancel_mode = 2;
numChirp = size(framesData,1);
%--------------------------------均值对消
if(cancelFlag == meanCancel_mode)
    meanCancelData = framesData - repmat(mean(framesData,1),[numChirp 1]);%均值对消 去零频
    %  repmat 函数将这个均值矩阵扩展为与 framesData 相同的维度
    Cancel_PulseCompre_Data = ifft(meanCancelData,rangeFFTNum,2);
    %脉冲压缩通过对数据进行傅里叶逆变换（IFFT）来实现
    %脉冲压缩 X = ifft(Y,n,dim) 返回沿维度 dim 的傅里叶逆变换。例如，如果 Y 是矩阵，则 ifft(Y,n,2) 返回每一行的 n 点逆变换。
end
%--------------------------------两脉冲对消
if(cancelFlag == twopCancel_mode)
    twopCancelData = framesData(1:end-pluseCancelSpace,:,:) - framesData(1+pluseCancelSpace:end,:,:);
    Cancel_PulseCompre_Data = ifft(twopCancelData,rangeFFTNum,2);%脉冲压缩
end
%--------------------------------三脉冲对消
if(cancelFlag == thrpCancel_mode)
    thrpCancelData = framesData(1:end-2,:,:) - 2 * framesData(2:end-1,:,:) + framesData(3:end,:,:);
    Cancel_PulseCompre_Data = ifft(thrpCancelData,rangeFFTNum,2);%脉冲压缩
end