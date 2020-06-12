function [PTEvM,pvalvM] = connectivity_basedPTEvariants(xM,variant,nc,h,m,tau,nnei,nsur)
% function [PTEvM,pvalvM] = connectivity_basedPTEvariants(xM,variant,nc,h,m,tau,nnei,nsur)
% Estimation of the selected connectivity-based PTE variants:
% '1A', '1B', '2A', '2B', '3A', '3B', '4A', '4B'
%INPUTS
%- xM     : multivariate time series (n x K)
%- variant: '1A' or '1B'or '2A' or '2B' or '3A' or '3B' or '4A' or '4B'
%- nc     : number of
%- h      : step ahead
%- m      : embedding dimension
%- tau    : lag
%- nnei   : number of neighbors
%- nsur   : number of surrogates
%OUTPUTS
%- PTEM   : PTE from selected variant for all posible directions, size K x K
%          (rows drive columns)
%- pvalM  : p-values from one-sided test using time shifted surrogates 
%           for all posible directions, size K x K
% Code written by Angeliki Papana (University of Macedonia, Greece)


if nargin ==6;
    nnei = 10;
    nsur = 100;
elseif nargin ==7;
    nsur = 100;
end

[~,K] = size(xM);

if floor(nc) ~= ceil(nc) || nc < 0 || nc > K
    error('nc must be a possitive integer smaller than K')
end

PTEvM = NaN*ones(K,K);  % estimated PTE variants for all pairs of directions
pvalvM = NaN*ones(K,K); % estimated p-values from PTE variants for all pairs of directions

for ik = 1:K
    for jk = ik+1:K
        % 1) xM(:,ik) -> xM(:,jk)| subset of xM(:,restV)

        % conditioning variables
        restV = setxor([ik jk]',(1:K)');
        
        if strcmp(variant,'1A')
            [zM1v,nlcv1] = subzM1A(xM(:,ik),xM(:,restV),nc);            
        elseif strcmp(variant,'1B')
            [zM1v,nlcv1] = subzM1B(xM(:,ik),xM(:,restV),nc);
        elseif strcmp(variant,'2A')
            [zM1v,nlcv1] = subzM2A(xM(:,ik),xM(:,restV),nc);            
        elseif strcmp(variant,'2B')
            [zM1v,nlcv1] = subzM2B(xM(:,ik),xM(:,restV),nc);
        elseif strcmp(variant,'3A')
            [zM1v,nlcv1] = subzM1A(xM(:,jk),xM(:,restV),nc);
        elseif strcmp(variant,'3B')
            [zM1v,nlcv1] = subzM1B(xM(:,jk),xM(:,restV),nc);
        elseif strcmp(variant,'4A')            
            [zM1v,nlcv1] = subzM2A(xM(:,jk),xM(:,restV),nc);
        elseif strcmp(variant,'4B')            
            [zM1v,nlcv1] = subzM2B(xM(:,jk),xM(:,restV),nc);
        end
        
        if nlcv1 == 0,  
            [PTEvM(ik,jk),surM1] = TExy(xM(:,ik),xM(:,jk),h,m,tau,nnei,nsur);
        else
            [PTEvM(ik,jk),surM1] = PTEXYZ(xM(:,ik),xM(:,jk),zM1v,h,m,tau,nnei,nsur);
        end
        
        % 2) xM(:,jk) -> xM(:,ik)| subset of xM(:,restV)
        if strcmp(variant,'1A')            
            [zM2v,nlcv2] = subzM1A(xM(:,jk),xM(:,restV), nc);
        elseif strcmp(variant,'1B')
            [zM2v,nlcv2] = subzM1B(xM(:,jk),xM(:,restV),nc);
        elseif strcmp(variant,'2A')
            [zM2v,nlcv2] = subzM2A(xM(:,jk),xM(:,restV),nc); 
        elseif strcmp(variant,'2B')
            [zM2v,nlcv2] = subzM2B(xM(:,jk),xM(:,restV),nc);            
        elseif strcmp(variant,'3A')
            [zM2v,nlcv2] = subzM1A(xM(:,ik),xM(:,restV),nc);
        elseif strcmp(variant,'3B')
            [zM2v,nlcv2] = subzM1B(xM(:,ik),xM(:,restV),nc);
        elseif strcmp(variant,'4A')
            [zM2v,nlcv2] = subzM2A(xM(:,ik),xM(:,restV),nc);
        elseif strcmp(variant,'4B')
            [zM2v,nlcv2] = subzM2B(xM(:,ik),xM(:,restV),nc);
        end
        
        if nlcv2 == 0,  
            [PTEvM(jk,ik),surM2] = TExy(xM(:,jk),xM(:,ik),h,m,tau,nnei,nsur);
        else
            [PTEvM(jk,ik),surM2] = PTEXYZ(xM(:,jk),xM(:,ik),zM2v,h,m,tau,nnei,nsur);
        end
        
        % estimate p-values
        pvalvM(ik,jk) = resampledonesidedpvalue([PTEvM(ik,jk);surM1]);
        pvalvM(jk,ik) = resampledonesidedpvalue([PTEvM(jk,ik);surM2]);
    end
end