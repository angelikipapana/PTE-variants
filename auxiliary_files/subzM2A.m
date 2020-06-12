function [zM2A,nlcv] = subzM2A(xV,zM,nc)
% function [zM2A,nlcv] = subzM2A(xV,zM)
% Estimate the subset of nc variables from zM, so that they are the least
% correlated (linearly) to xV (based on the criteria 2A, estimation based
% on correlation coefficient)
% if xV is correlated with variables from zM, the subset is []
%INPUTS
% - xV : driving variable
% - zM : set of (conditioning) variables
% - nc : number of conditioning variables that should be chosen from zM
%OUTPUTS
% zM2A : subset of zM with nc variables in total (criteria 2A)
% nlcv: number of variables non-significantly (linearly) correlated with the
%       driving one


% Compute sample correlation & p-values between
% xV (driving variable) and zM (conditioning variables)
[rM,pM] = corrcoef([xV zM]);

% find least significant correlations with driving variable and
% find a reduced number of conditioning variables:
% I consider only the 1st line (associations of 1rst variable (driving)
% with all the other, deleate the 1rst element of diagonal at position 1x1:
crM = abs(rM(1,2:end)); % estimated correlation
cpM = pM(1,2:end);      % estimated p-values

% Find non-significant p-values (>=0.05)
[~,sigb] = find(cpM>=0.05);

if isempty(sigb)
    zM2A = [];
    nlcv = 0;
else
    
    % Only non-significant p-values
    sigrM = crM(1,sigb);
    
    % Sort in ascending order the estimated correlation with non-signif p-values
    [~,bx] = sort(sigrM,'ascend');
    
    % bx: the 'order' of the non-significant correlations (ascending order)
    cV = sigb(bx);
    
    % number of variables non-significantly correlated with the driving one
    nlcv = length(cV);
    
    if nlcv >= nc,
        zM2A = zM(:,cV(1:nc));
    elseif nlcv < nc && nlcv ~=0
        zM2A = zM(:,cV);
    end
end
