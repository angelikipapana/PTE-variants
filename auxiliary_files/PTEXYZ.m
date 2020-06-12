function [pTExyz,surPTEx] = PTEXYZ(xV,yV,zM,h,m,tau,nnei,nsur)
% [PTEXY,surM] = PTExy(xV,yV,zM,h,m,tau,nnei,nsur)
% This function estimates the Partial Transfer Entropy (PTE) 
% based on the k-nearest neighbors (KNN) estimator (Kraskov et al., 2004)
% Also, time shifted surrogates are generated for the driving variable and the
% corresponding surrogate PTE values are exported
%INPUTS
% - xV    : driving time series, size: n x 1
% - yV    : response time series, size: n x 1
% - zM    : conditioning time series, size: n x s
% - h     : h steps ahead 
% - m     : embenddings dimension  
% - tau   : lag
% - nnei  : number of neighbors
%OUTPUTS
% - pTExyz  : partial transfer entropy, PTE(X->Y|Z)
% - surPTEx : surrogate PTE value, PTE(X*->Y|Z)
%             X*: time-shifted time series of X  
% Code written by Angeliki Papana (University of Macedonia, Greece)  

% rescale time series to be in [0,1] 
xV = rangescale(xV); 
yV = rangescale(yV);
zM = rangescale(zM);

[N,s] = size(zM); 
M = (m-1)*tau;

yVh = NaN*ones(N-h-M,1);

xM = NaN*ones(N-h-M,m);
yM = NaN*ones(N-h-M,m);
dzM = NaN*ones(N-h-M,s*m);
 
for i = M+1:N-h 
    yVh(i-M,1) = yV(i+h);
    
    xM(i-M,:) = xV(i-(m-1)*tau:tau:i);
    yM(i-M,:) = yV(i-(m-1)*tau:tau:i);
    for is = 1:s
        dzM(i-M,(is-1)*m+1:is*m) = zM(i-(m-1)*tau:tau:i,is);
    end
end

pTExyz=cmikra1n(yVh,xM,[yM dzM],nnei);   % X->Y|Z

% Surrogate PTE values - time shifted surrogates
offsetfrac = 0.1; % For time-shifted surrogates
offset = round(N*offsetfrac/2); % time offset for surrogates

surPTEx = NaN*ones(nsur,1);
for isur = 1:nsur
    istart = offset+unidrnd(N-2*offset);
    sxV = xV([istart:N 1:istart-1],:);
    sxM = NaN*ones(N-h-M,m);   
    for i = M+1:N-h 
        sxM(i-M,:) = sxV(i-(m-1)*tau:tau:i);       
    end   
    surPTEx(isur,1)=cmikra1n(yVh,sxM,[yM dzM],nnei);   % X*->Y|Z
end