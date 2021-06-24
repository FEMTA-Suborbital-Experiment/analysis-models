all_data=importdata('cells.out');
x=all_data(:,1);
y=all_data(:,2);
z=all_data(:,3);

scatter3(x,y,z,'.');