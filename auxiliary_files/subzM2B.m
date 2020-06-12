function [zM2B,nlcv] = subzM2B(xV,zM,nc,nnei,nsur)
% function [zM2B,nlcv] = subzM2B(xV,zM)
% Estimate the subset of nc variables from zM, so that they are the least
% correlated (nonlinearly) to xV based on the criteria 2B, i.e.  
% in terms of mutual information
% if xV is correlated with ALL variables from zM, the subset is [].
%INPUTS
% - xV : driving variable 
% - zM : set of (conditioning) variables 
% - nc : number of conditioning variables that should be chosen from zM
%OUTPUTS
% zM2B : subset of zM with nc variables in total (criteria 2B)
% nlcv: number of variables not-significantly (nonlinearly) correlated with the 
%       driving one

if nargin ==3;
    nnei = 10;
    nsur = 100;
end

% subset of zM based on condition 2B
% zM2B = NaN*ones(length(xV),nc);

%ncv: number of variables in zM
[n,ncv]=size(zM);  

% Compute mutual information & p-values between 
% xV (driving variable) and zM (conditioning variables)
rM = NaN*ones(1,ncv);
for in=1:ncv
    rM(1,in) = mkraskov1(xV,zM(:,in),nnei);    
end
crM = abs(rM);  % estimated absolute correlation


% Random permutation of driving time series and estimation of 'surrogate' 
% MI for nsur 'surrogates'
rrM = NaN*ones(nsur,ncv);
for isur = 1:nsur
    rpV = randperm(n);
    rxV = xV(rpV); 
    for jn=1:ncv
        rrM(isur,jn) = mkraskov1(rxV,zM(:,jn),nnei);
    end
end

% p-values from one sided-test 
cpM = NaN*ones(1,ncv);
for ik=1:ncv
    cpM(1,ik) = resampledonesidedpvalue([rM(1,ik); rrM(:,ik)]);
end

% Find non-significant p-values (>=0.05)
[~,sigb] = find(cpM>=0.05);

if isempty(sigb)
    zM2B = [];
    nlcv = 0;
else
    % Only estimated correlation with non-significant p-values
    sigrM = crM(1,sigb);
    
    % Sort in ascending order the estimated correlation with non-signif p-values
    [~,bx] = sort(sigrM,'ascend');
    
    % bx: the 'order' of the non-significant correlations (ascending order)
    cV = sigb(bx);
    
    % number of variables significantly correlated with the driving one
    nlcv = length(cV);
    
    if nlcv >= nc,
        zM2B = zM(:,cV(1:nc));
    elseif nlcv < nc && nlcv ~= 0
        zM2B = zM(:,cV);
    end
end
