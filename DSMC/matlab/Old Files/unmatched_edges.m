alldata=load('unmatched edges.txt');
x1=alldata(:,1);
scatter3(alldata(:,1),alldata(:,2),alldata(:,3),'r.')
hold on
scatter3(alldata(:,4),alldata(:,5),alldata(:,6),'r.')
xlabel('x')
ylabel('y')
zlabel('z')

all_data=load('pointsMatlab.txt');
x=all_data(:,2);
y=all_data(:,3);
z=all_data(:,4);

% 
% indicesx=find(x<250);
% x(indicesx)=0;
% indicesy=find(y<40);
% y(indicesy)=0;
% indicesz=find(z>10);
% z(indicesz)=0;

scatter3(x,y,z,'bo')
xlabel('x')
ylabel('y')
zlabel('z')
