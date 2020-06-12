function yM=rangescale(xM)
n1=size(xM,1);
yM=(xM-repmat(min(xM),n1,1))./repmat(range(xM),n1,1);