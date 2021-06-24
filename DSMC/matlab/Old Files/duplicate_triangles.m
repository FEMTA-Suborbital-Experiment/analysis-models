x=load('triangles.surf');
[u,I,J]=unique(x,'rows','first');
hasDuplicates=size(u,1) < size(x,1);
ixDupRows=setdiff(1:size(x,1),I);
dupRowValues=x(ixDupRows,:);