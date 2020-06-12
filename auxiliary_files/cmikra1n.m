function [mi]=cmikra1n(xM,yM,zM,k)
% function [mi]=cmikra1n(xM,yM,zM,k)
% Estimated Conditional Mutual Information mi = I(xM,yM|zM)
 

n1=size(xM,1);
psimar=NaN*ones(n1,3);
xM=rangescale(xM);
yM=rangescale(yM);
zM=rangescale(zM);
xMb=[xM yM zM];

% for 64-bit windows
[nnidx, dists] = annMaxquery(xMb', xMb', k+1);

maxdistV=dists(end,:)';

nz=npoinmultranges(zM,maxdistV-ones(n1,1)*10^(-10));
nyz=npoinmultranges([yM zM],maxdistV-ones(n1,1)*10^(-10));
nxz=npoinmultranges([xM zM],maxdistV-ones(n1,1)*10^(-10));

psimar(:,1)=psi(nxz);
psimar(:,2)=psi(nyz);
psimar(:,3)=-1*psi(nz);

mi=psi(k)-mean(sum(psimar,2));

