clear
clc
alldata=load('femta_matlab.surf');
xmin=min(alldata(:,2))/500;
xmax=max(alldata(:,2))/500;
ymin=min(alldata(:,3))/500;
ymax=max(alldata(:,3))/500;
zmin=min(alldata(:,4))/500;
zmax=max(alldata(:,4))/500;

[xmin xmax; ymin ymax; zmin zmax]