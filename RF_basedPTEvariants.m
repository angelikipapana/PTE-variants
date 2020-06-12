function [PTEvM,pvalvM,result] = RF_basedPTEvariants(variant,datname,path,datdir,savedir,savename1,savename2,h,m,tau,nnei,nsur)
%function [PTEvM,pvalvM,result] = RF_basedPTEvariants(variant,datname,path,datdir,savedir,savename1,savename2,h,m,tau,nnei,nsur)
% Estimation of the selected RF-based PTE variant from '5A', '5B', '5C' 
%INPUTS
%- variant     : '5A' or '5B' or '5C' multivariate time series (n x K)
%- datname     : file with data, a multivariate time series of size n (length) x K (number
%                of variables), e.g. datname = 'test.dat'; 
% - path       : directory of the Rscript 'optimalSubsetGenerator.R'
%               e.g. path = 'C:/PTEvariants/';
% - datdir     : set the directory of data (file datname)
% - savedir    : set the directory for saving the outputs with the depths
%                and ranks from RF 
% - savename1  : file with ranks from RF, e.g. 'ranks.dat' 
% - savename2  : file with depths from RF, e.g. 'depths.dat' 
% - h      : step ahead
% - m      : embedding dimension
% - tau    : lag
% - nnei   : number of neighbors
% - nsur   : number of surrogates
%OUTPUTS
%- PTEM   : PTE from selected variant for all posible directions, size K x K
%          (rows drive columns)
%- pvalM  : p-values from one-sided test using time shifted surrogates 
%           for all posible directions, size K x K
% Codes written by Angeliki Papana (University of Macedonia, Greece) &
% Ariadni Papana-Dagiasis (Cleveland State University, USA)
% Created: 1/6/2020

if nargin ==10;
     nnei = 10;     % number of neighbors
     nsur = 100;    % number of surrogates
elseif nargin ==11;
     nsur = 100;
end

xM = importdata(strcat(datdir,datname));
[n,K] = size(xM);

PTEvM = NaN*ones(K,K);  % estimated PTE variants for all pairs of directions
pvalvM = NaN*ones(K,K); % estimated p-values from PTE variants for all pairs of directions

if exist(strcat(savedir,savename1), 'file')==2
  delete(strcat(savedir,savename1));
end
if exist(strcat(savedir,savename2), 'file')==2
  delete(strcat(savedir,savename2));
end


if strcmp(variant,'5A') || strcmp(variant,'5B') 
    system(strcat(['Rscript ',path,'optimalSubsetGenerator.R ',datdir,datname,' ',savedir,savename1,' ',savedir,savename2]));
elseif strcmp(variant,'5C')
   system(strcat(['Rscript ',path,'optimalSubsetGenerator5C.R ',datdir,datname,' ',savedir,savename1,' ',savedir,savename2]));
end

result = importdata(strcat([savedir,savename1]));

PTEvM = NaN*ones(K,K);  % PTE from case
pvalvM = NaN*ones(K,K); % p-values from case

for ik = 1:K
    for jk = ik+1:K        
         % 1) xM(:,ik) -> xM(:,jk)| subset of xM from 5A or 5B or 5C
        if strcmp(variant,'5A')
            vect1 = result(ik,2:end);
            vect2 = setdiff(vect1,0);
            fcondvar = setdiff(vect2,jk);
        elseif strcmp(variant,'5B')
            vect1 = result(jk,2:end);
            vect2 = setdiff(vect1,0);
            fcondvar = setdiff(vect2,ik);
        elseif strcmp(variant,'5C')
            a5c = find(result(:,1)==ik & result(:,2)==jk);
            vect1 = result(a5c,3:end);
            fcondvar = setdiff(vect1,0);
        end
        if isempty(fcondvar)
            [PTEvM(ik,jk),surM1] = TExy(xM(:,ik),xM(:,jk),h,m,tau,nnei,nsur);
        else
            zM1v = xM(:,fcondvar);
            [PTEvM(ik,jk),surM1] = PTEXYZ(xM(:,ik),xM(:,jk),zM1v,h,m,tau,nnei,nsur);
        end    

        % 2) xM(:,jk) -> xM(:,ik)| subset of xM from 5A or 5B or 5C
        if strcmp(variant,'5A')
            avect1 = result(jk,2:end);
            avect2 = setdiff(avect1,0);
            afcondvar = setdiff(avect2,ik);
        elseif strcmp(variant,'5B')
            avect1 = result(ik,2:end);
            avect2 = setdiff(avect1,0);
            afcondvar = setdiff(avect2,jk);
        elseif strcmp(variant,'5C')
            av = find(result(:,1)==jk & result(:,2)==ik);
            avect1 = result(av,3:end);
            afcondvar = setdiff(avect1,0);
        end
        if isempty(afcondvar)
            [PTEvM(jk,ik),surM2] = TExy(xM(:,jk),xM(:,ik),h,m,tau,nnei,nsur);
        else
            zM2v = xM(:,afcondvar);
            [PTEvM(jk,ik),surM2] = PTEXYZ(xM(:,jk),xM(:,ik),zM2v,h,m,tau,nnei,nsur);
        end
        % estimate p-values 
        pvalvM(ik,jk) = resampledonesidedpvalue([PTEvM(ik,jk);surM1]);
        pvalvM(jk,ik) = resampledonesidedpvalue([PTEvM(jk,ik);surM2]);
    end
end