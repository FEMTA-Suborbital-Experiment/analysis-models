all_data=load('cathodePressTemp_timestep800.out');
pressure=all_data(:,4);
torr=pressure./133;
d=1.25;
P_by_d=d*mean(torr);