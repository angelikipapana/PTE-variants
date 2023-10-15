function [mi topsi1 topsi2]=mkraskov1(xM1,xM2,k)
n1=size(xM1,1);


psis=psi(k)+psi(n1);



xMb=[xM1 xM2];

% [nnidx, dists] = annquerymax(xMb', xMb', k+1);
[nnidx, dists] = annMaxquery(xMb', xMb', k+1);
maxdistV=dists(end,:)';



topsi1=npoinmultranges(xM1,maxdistV-ones(n1,1)*10^(-10));


topsi2=npoinmultranges(xM2,maxdistV-ones(n1,1)*10^(-10));


psiff=psi([topsi1 topsi2]);

% %test
% holo=find(isfinite(sum(psiff,2)));
% psiff=psiff(holo,:);
% psiff=psi([topsi1(holo) topsi2(holo)]);
% %test

mi=psis-mean(sum(psiff,2));