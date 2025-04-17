function [retVal] = CFAR_2D(RDM)

    %% CFAR检测算法
    %通过完整的距离多普勒图滑动窗口
    %在两个维度中选择参考单元的数量
    Tr = 8; 
    Td = 4;
    %选择被测单元（CUT）周围两个维度的保护单元数量，以进行准确CFAR检测
    Gr = 4;
    Gd = 2;
    
    %原矩阵行列数
    %[r,c]=size(RDM);%原矩阵行列数
    [r,c]=size(RDM);
    CFAR = zeros(r,c);       %所以cfaf返回的是行列相同的结果
    
    %为参考单元上的每次迭代创建一个向量来存储noise_level（噪声电平）
    r_margin = 2*(Tr+Gr);
    d_margin = 2*(Td+Gd);
    
    
    gridSize = (2*Tr+2*Gr+1)*(2*Td+2*Gd+1);
    numGcells = (2*Gr+1)*(2*Gd+1);
    numTcells = gridSize - numGcells;           %参考单元数量
    
    
    pfa=1e-4;                                   %虚警概率  %需要注意的是这里的参数被我改了，原来是pfa=1e-6;
    alpha=numTcells*(pfa^(-1/numTcells)-1);     %门限系数
    snr_offset =4.5*log10(alpha);               %以 dB 为单位的 SNR 值偏移阈值
    
    RDM_temp1=zeros(r,c);
    for i=1:r
        for j=1:c
            RDM_temp1(i,j)=10*log10(abs(RDM(i,j)));
        end
    end
    
    a=mean(RDM_temp1,'all');
    RDM_temp=padarray(RDM_temp1,[Tr+Gr,Td+Gd],a,'both');          %采用补平均数法
    noise_level = zeros(r,c);
    
    for i = 1:r                                                   % 距离边界
        for j = 1:c                                               % 多普勒边界    
            sig_pow = (RDM_temp(i:i+r_margin,j:j+d_margin ));     %从db变为线性
            G_pow   = (RDM_temp(i+Tr:i+Tr+Gr*2,j+Td:j+Td+Gd*2));  %保护单元
            sig_sum = sum(sum(sig_pow))-sum(sum(G_pow));          % 所有参考单元内的总和信号
            noise_level(i,j) = (sig_sum/numTcells);               % 对所有使用的参考单元格的总和值求平均
            sig_threshold = noise_level(i,j) + snr_offset;        % 添加偏移量
            
            % 将 CUT 与阈值进行比较
            if (RDM_temp(i+r_margin/2, j+d_margin/2) > sig_threshold)
                CFAR(i,j) = RDM_temp1(i,j);
            end  
        end
    end
    
    retVal=CFAR;
end