all_data=load('points.txt');
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

scatter3(x,y,z)
xlabel('x')
ylabel('y')
zlabel('z')