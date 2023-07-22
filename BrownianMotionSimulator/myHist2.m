function [binMids, NormCounts] = myHist2(data, binW)
edge = -1:binW:1;
[NormCounts,edges] = histcounts(data,edge,'Normalization','probability');
binMids = edges(2:end) - (edges(2)-edges(1))/2;
end
