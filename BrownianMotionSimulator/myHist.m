function [binMids, NormCounts, totalCounts] = myHist(data, binW)
%This is a function that will return the histogram counts and the midpoints
%of the histogram bins in order to plot a "smooth" histogram.
%Input: data- single vector of data to make the histogram
%       binW- this is how wide to make each bin
%Output: binMids- vector of the midpoints of the bins
%        NormCounts- vector of normalized portion of how many counts (from 
%                    the data) are in each bin portion
%        TotCounts- this is the total number of counts from data that is
%                   used to do the normalization.

h = histogram(data, 'BinWidth', binW);

counts = h.Values;
totalCounts = sum(counts);
NormCounts = counts/totalCounts;
binEdges = h.BinEdges;

binMids = movmean(binEdges, 2);
binMids = binMids(2:end);

end