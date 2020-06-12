function [PTEM,surM] = PTE(xV,yV,zM,h,m,tau,k,nsur)
% function [PTEM,surM] = PTE(xV,yV,zM,h,m,tau,k,nsur)
% This function estimates the Partial Transfer Entropy (PTE) 
% based on the k-nearest neighbors (KNN) estimator (Kraskov et al., 2004)
% Also, time shifted surrogates are generated for the driving variable and the
% corresponding surrogate PTE values are exported
%INPUTS
% - xV,yV : time series, each of size n x 1
% - zM    : s time series (size: n x s) 
% - h     : h steps ahead 
% - m     : embedding dimension  
% - tau   : lag
% - k     : number of neighbors
%OUTPUTS
% - PTEM : partial transfer entropy PTEM = [PTE(X->Y|Z) PTE(Y->X|Z)] (size 2 x 1)
% - surM : surrogate PTE surM = [PTE(X*->Y|Z) PTE(Y*->X|Z)]
%          X*, Y*: time-shifted time series of X and Y, respectively 
% Code written by Angeliki Papana (University of Macedonia, Greece)  

% rescale time series to be in [0,1] 
xV = rangescale(xV); 
yV = rangescale(yV);
zM = rangescale(zM);

[N,s] = size(zM); %N = length(xV);
M = (m-1)*tau;

xVh = NaN*ones(N-h-M,1);
yVh = NaN*ones(N-h-M,1);

xM = NaN*ones(N-h-M,m);
yM = NaN*ones(N-h-M,m);
dzM = NaN*ones(N-h-M,s*m);
 
for i = M+1:N-h 
    xVh(i-M,1) = xV(i+h);
    yVh(i-M,1) = yV(i+h);
    
    xM(i-M,:) = xV(i-(m-1)*tau:tau:i);
    yM(i-M,:) = yV(i-(m-1)*tau:tau:i);
    for is = 1:s
        dzM(i-M,(is-1)*m+1:is*m) = zM(i-(m-1)*tau:tau:i,is);
    end
end

pTExyz=cmikra1n(yVh,xM,[yM dzM],k);   % X->Y|Z
pTEyxz=cmikra1n(xVh,yM,[xM dzM],k);   % Y->X|Z

PTEM = [pTExyz pTEyxz];

% Surrogate PTE values - time shifted surrogates
offsetfrac = 0.1; % For time-shifted surrogates
offset = round(N*offsetfrac/2); % time offset for surrogates

surPTEx = NaN*ones(nsur,1);
surPTEy = NaN*ones(nsur,1);
for isur = 1:nsur
    istart = offset+unidrnd(N-2*offset);
    sxV = xV([istart:N 1:istart-1],:);
    syV = yV([istart:N 1:istart-1],:);
 
    sxM = NaN*ones(N-h-M,m);
    syM = NaN*ones(N-h-M,m);
   
    for i = M+1:N-h 
        sxM(i-M,:) = sxV(i-(m-1)*tau:tau:i);       
        syM(i-M,:) = syV(i-(m-1)*tau:tau:i);       
    end
    
    surPTEx(isur,1)=cmikra1n(yVh,sxM,[yM dzM],k);   % X->Y|Z
    surPTEy(isur,1)=cmikra1n(xVh,syM,[xM dzM],k);   % Y->X|Z
end

surM = [surPTEx surPTEy];
