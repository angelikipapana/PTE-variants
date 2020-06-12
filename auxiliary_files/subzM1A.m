function [zM1A,nlcv] = subzM1A(xV,zM,nc)
% function [zM1A,nlcv] = subzM1A(xV,zM)
% Estimate the subset of nc variables from zM, so that they are the most
% correlated (linearly) to xV
% for PTE variant 1A, estimation based on correlation coefficient
% if xV is not correlated with variables from zM, the subset is []
%INPUTS
% - xV  : driving variable
% - zM  : set of (conditioning) variables
% - nc  : number of conditioning variables that should be chosen from zM
%OUTPUTS
% -zM1A : subset of zM with nc variables in total (criteria 1A)
% -nlcv : number of variables significantly (linearly) correlated with the
%       driving one
% Code written by Angeliki Papana (University of Macedonia, Greece)


% Compute sample correlation & p-values between
% xV (driving variable) and zM (conditioning variables)
[rM,pM] = corrcoef([xV zM]);

% find most significant correlations with driving variable
%To find a reduced number of conditioning variables:
% I consider only the 1st line (associations of 1rst variable (driving)
% with all the other, deleate the 1rst element of diagonal at position 1x1:
crM = abs(rM(1,2:end)); % estimated correlation
cpM = pM(1,2:end);      % estimated p-values

% Find significant p-values (<0.05)
[~,sigb] = find(cpM<0.05);

if isempty(sigb)
    zM1A = [];
    nlcv=0;	
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
        zM1A = zM(:,cV(1:nc));
    elseif nlcv < nc && nlcv ~= 0
        zM1A = zM(:,cV);
    end
end
