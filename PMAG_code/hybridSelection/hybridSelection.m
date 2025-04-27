%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the code for the U-SPEC algorithm, which is proposed in   %
% the following paper:                                              %
%                                                                   %
% D. Huang, C.-D. Wang, J.-S. Wu, J.-H. Lai, and C.-K. Kwoh.        %
% "Ultra-Scalable Spectral Clustering and Ensemble Clustering."     %
% IEEE Transactions on Knowledge and Data Engineering, 2020.        %
% DOI: https://doi.org/10.1109/TKDE.2019.2903410                    %
%                                                                   %
% The code has been tested in Matlab R2016a and Matlab R2016b.      %
% Website: https://www.researchgate.net/publication/330760669       %
% Written by Huang Dong. (huangdonghere@gmail.com)                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RpFea = hybridSelection(fea, p, distance)

if nargin < 3
    distance = 'euclidean'; % Use Euclidean distance by default.
end

N = size(fea,1);
if p>N
    p = N;
end

warning('off');

%% Get $p$ representatives by hybrid selection
RpFea = getRepresentativesByHybridSelection(fea, p, distance);

function RpFea = getRepresentativesByHybridSelection(fea, pSize, distance, cntTimes)
% Huang Dong. Mar. 20, 2019.
% Select $pSize$ representatives by hybrid selection.
% First, randomly select $pSize * cntTimes$ candidate representatives.
% Then, partition the candidates into $pSize$ clusters by k-means, and get
% the $pSize$ cluster centers as the final representatives.

if nargin < 4
    cntTimes = 10;
end

N = size(fea,1);
bigPSize = cntTimes*pSize;
if pSize>N
    pSize = N;
end
if bigPSize>N
    bigPSize = N;
end

rand('state',sum(100*clock)*rand(1));
bigRpFea = getRepresentivesByRandomSelection(fea, bigPSize);

if strcmp(distance,'euclidean')
    [~, RpFea] = litekmeans(bigRpFea,pSize,'MaxIter',10);
else
    [~, RpFea] = litekmeans(bigRpFea,pSize,'MaxIter',10,'Distance',distance);
end

function [RpFea,selectIdxs] = getRepresentivesByRandomSelection(fea, pSize)
% Huang Dong. Mar. 20, 2019.
% Randomly select pSize rows from fea.

N = size(fea,1);
if pSize>N
    pSize = N;
end
selectIdxs = randperm(N,pSize);
RpFea = fea(selectIdxs,:);
