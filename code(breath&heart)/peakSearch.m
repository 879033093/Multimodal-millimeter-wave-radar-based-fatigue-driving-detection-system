function [RD_pearkSearch,peakSearchList] = peakSearch(RD_cfar,cfarTargetList)
%注意cfarTargetList中第一行包含距离索引，第二行包含相应的多普勒索引

    peakSearchList = [];
    [rangeLen, dopplerLen] = size(RD_cfar);
    RD_pearkSearch = zeros(rangeLen, dopplerLen);
    length = size(cfarTargetList,2); %(cfarTargetList,2)表示cfarTargetList的列，值得注意的是：cfarTargetList的列恰好是每一个检测点的坐标
    %因此length代表了一共检测到了多少个目标点
    for targetIdx = 1:length   %对所有检测点进行遍历

        rangeIdx   = cfarTargetList(1,targetIdx); %提取出每一个点的距离坐标
        dopplerIdx = cfarTargetList(2,targetIdx); %提取出每一个点的速度坐标

        if rangeIdx > 1 && rangeIdx < rangeLen && dopplerIdx > 1 && dopplerIdx < dopplerLen %边界点不考虑  
           
            if RD_cfar(rangeIdx,dopplerIdx) > RD_cfar(rangeIdx - 1,dopplerIdx) && ...
                    RD_cfar(rangeIdx,dopplerIdx) > RD_cfar(rangeIdx + 1,dopplerIdx) && ...
                    RD_cfar(rangeIdx,dopplerIdx) > RD_cfar(rangeIdx,dopplerIdx - 1) && ...
                    RD_cfar(rangeIdx,dopplerIdx) > RD_cfar(rangeIdx,dopplerIdx + 1)         %这四行就是在表示找一个点，他比上下左右都大

                    RD_pearkSearch(rangeIdx,dopplerIdx) = RD_cfar(rangeIdx,dopplerIdx);% RD_pearkSearch为RD_cfar在（rangeIdx,dopplerIdx)下的值

                    cfarTarget = [rangeIdx ; dopplerIdx];

                    peakSearchList = [peakSearchList cfarTarget];
            end   
        end
    end  
end
