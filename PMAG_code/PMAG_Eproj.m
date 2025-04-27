function [W,alpha,obj,runtime] = PMAG_Eproj(X,projDim,c,r)
% This function is a modification of the code provided by the following work.
% Ref:
% Qianyao Qiang, Bin Zhang, Jason Chen Zhang, Chaodie Liu, Feiping Nie.
% Projection with Mixed-Size Anchor Graphs.
% X: n*d data matrix
% projDim: project dimension
tic;
NITER = 30;
X =(double(X));%n*d
X = mapstd(X',0,1);
% projDim = c;
%% Initalization
numNearestAnchor = 5;
anchorGraphNum = 4;
B = cell(1,anchorGraphNum);
for index = 1:anchorGraphNum
    index1 = index*10;%*3
    numAnchor = index1*c;
    locAnchor = hybridSelection(X', numAnchor)'; locAnchor = locAnchor';% hybridSelection
    Z = ConstructA_NP(X,locAnchor',numNearestAnchor);
    sumZ = sum(Z);
    sqrtZ = sumZ.^(-0.5);
    B{1,index} = Z*(diag(sqrtZ));
end

alpha = ones(anchorGraphNum,1)/anchorGraphNum;

% begin the loop
obj = [];
iter = 1;
St = X*X';
while iter < NITER
    %% update W
    B1 = [];
    for a = 1:anchorGraphNum
        B1 = cat(2, B1, sqrt(alpha(a)^r) * B{1,a});
    end
    XB = X*B1;
    XB2 = XB*XB';    
    alphaSum = sum(alpha.^r);
    Wm = (alphaSum*St+eps)\XB2;
    [W,~,~] = eig1(Wm,projDim,1,0);
    
    %% update alpha
    gs = zeros(anchorGraphNum,1);
    WX = W'*X;
    WX2 = WX*WX';
    for a = 1:anchorGraphNum
        WXB = WX*B{1,a};
        WXB2 = WXB*WXB';
        gs(a) = trace(WX2-WXB2);
    end
    gsw = gs;
    if r > 1
        alpha = gs.^(1/(1-r));
        alpha = alpha./sum(alpha);
    else
        error('r must larger than 1.\n');
    end
    
    %% calculate objective value
    obj(1,iter) = sum((alpha.^r).*gsw);    
    if iter>2 && abs( obj(iter-1)-obj(iter)) < 1e-10
        break;
    end
    iter = iter + 1;
end
runtime = toc;
end

