function [zM1B,nlcv] = subzM1B(xV,zM,nc,nnei,nsur)
% function [zM1B,nlcv] = subzM1B(xV,zM)
% Estimate the subset of nc variables from zM, so that they are the most
% correlated (nonlinearly) to xV
% for PTE variant 1B, i.e. in terms of mutual information
% if xV is not correlated with variables from zM, the subset is [].
%INPUTS
% - xV  : driving variable
% - zM  : set of (conditioning) variables
% - nc  : number of conditioning variables that should be chosen from zM
%OUTPUTS
% -zM1B : subset of zM with nc variables in total (criteria 1B)
% -nlcv : number of variables significantly (nonlinearly) correlated with the
%       driving one
% Code written by Angeliki Papana (University of Macedonia, Greece)


if nargin ==3;
    nnei = 10;
    nsur = 100;
end

%ncv: number of variables in zM
[n,ncv]=size(zM);

% Compute mutual information & p-values between
% xV (driving variable) and zM (conditioning variables)
rM = NaN*ones(1,ncv);
for in=1:ncv
    rM(1,in) = mkraskov1(xV,zM(:,in),nnei);
end

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

% find most significant correlations with driving variable
crM = abs(rM);  %: estimated absolute correlation

% Find significant p-values (<0.05)
[~,sigb] = find(cpM<0.05);

if isempty(sigb)
    zM1B = [];
    nlcv = 0;
else
    % Only estimated correlation with significant p-values
    sigrM = crM(1,sigb);
    
    % Sort in descending order the estimated correlation with signif p-values
    [~,bx] = sort(sigrM,'descend');
    
    % bx: the 'order' of the significant correlations (descending order)
    cV = sigb(bx);
    
    % number of variables significantly correlated with the driving one
    nlcv = length(cV);
    
    if nlcv >= nc,
        zM1B = zM(:,cV(1:nc));
    elseif nlcv < nc  && nlcv ~=0
        zM1B = zM(:,cV);
    end
end