function npV=npoinmultranges(xM,rV)

% for 64-bit windows
npV = annMaxRvaryquery(xM', xM',rV, 1, 'search_sch', 'fr', 'radius', sqrt(1));


npV=double(npV);
npV(npV==0)=1;


