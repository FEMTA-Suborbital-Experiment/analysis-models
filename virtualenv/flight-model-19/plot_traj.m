clear all;
close all;
clc;

load altitude
load accel
load time

hold on;
grid on;
yyaxis left
ytickformat('%,i');
ylabel('Altitude [km]');
plot(t ./ 60, h ./ 1000);
yyaxis right;
ytickformat('%i');
ylabel('Acceleration [g]');
plot(t ./ 60, a ./ 9.81);
xlabel('Time [min]');
ax = gca;
ax.FontSize = 14; 
