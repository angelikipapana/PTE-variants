function [TEx,sTEx] = TExy(xV,yV,h,m,tau,nnei,nsur)
% function [TEx,sTEx] = TExy(xV,yV,h,m,tau,nnei,nsur)
% This function estimates the bivariate Transfer Entropy (TE) from the 
% observed time series TE(X->Y) 
% and the TE from time shifted surrogates of the driving variable X
% TE(X*->Y), where X* is a time shifted surrogate time series from X 
% The k-nearest neighbor (KNN) estimator is used (Kraskov et al., 2004)
%INPUTS
% - xV            : driving time series, size: n x 1
% - yV            : response time series, size: n x 1
% - h             : h steps ahead of xV
% - m             : embenddings dimension 
% - tau           : lag 
% - nnei          : number of neighbors
% - nsur          : number of surrogates
% OUTPUTS
% - TExy          : TE(X->Y), size: 1x1
% - sTEx          : TE(X*->Y), size: nsur x 1  
%                   X*: time-shifted time series of X  
%
% Code written by Angeliki Papana (University of Macedonia, Greece)  

% rescale time series to be in [0,1] 
xV = rangescale(xV); 
yV = rangescale(yV);

N = length(xV);
M = (m-1)*tau;

yVh = NaN*ones(N-h-M,1);
xM = NaN*ones(N-h-M,m);
yM = NaN*ones(N-h-M,m);

for i = M+1:N-h 
    yVh(i-M,1) = yV(i+h);
    xM(i-M,:) = xV(i-(m-1)*tau:tau:i);
    yM(i-M,:) = yV(i-(m-1)*tau:tau:i);
end

TEx = cmikra1n(yVh,xM,yM,nnei);          % TE(X->Y)


% Surrogate PTE values - time shifted surrogates
offsetfrac = 0.1; % For time-shifted surrogates
offset = round(N*offsetfrac/2); % time offset for surrogates

sTEx  = NaN*ones(nsur,1);

for isur = 1:nsur
 istart = offset+unidrnd(N-2*offset);
    sxV = xV([istart:N 1:istart-1],:);
    sxM = NaN*ones(N-h-M,m);
    for i = M+1:N-h 
        sxM(i-M,:) = sxV(i-(m-1)*tau:tau:i);       
    end
  
  sTEx(isur,1) = cmikra1n(yVh,sxM,yM,nnei);  % sur TEx values
end
