function [pval] = resampledonesidedpvalue(xV)
% function [pval] = resampledonesidedpvalue(xV)
% Estimates p-values for one sided test
% xV = [original value; surrogate values]
% correction of Yu & Huang (2001)

nx = length(xV);
[oxV,ixV]=sort(xV);
rnkx = find(ixV == 1);
isamexV = find(xV==xV(1));
if length(isamexV)==nx
    rnkx=round(nx/2);
elseif length(isamexV)>=2
     irand = unidrnd(length(isamexV));
     rnkx = rnkx+irand-1;
end  
pval = (1-(rnkx-0.326)/(nx+0.348));
