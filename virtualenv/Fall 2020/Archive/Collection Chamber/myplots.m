plot((0:dt:totalTime),T290/1000,'b-','Linewidth',3)
hold on
plot((0:dt:totalTime),Pv290/1000,'b--','Linewidth',1)
plot((0:dt:totalTime),T300/1000,'g-','Linewidth',3)
plot((0:dt:totalTime),Pv300/1000,'g--','Linewidth',1)
plot((0:dt:totalTime),T310/1000,'r-','Linewidth',3)
plot((0:dt:totalTime),Pv310/1000,'r--','Linewidth',1)
ax=gca;
ax.FontSize=18;
title("Collection Chamber Pressure for Duration of Experiment",'Fontsize',22)
ylabel("Pressure [kPa]",'Fontsize',20)
xlabel("Time [s]",'Fontsize',20)
legend('Pressure at 290K','Vapor Pressure at 290K','Pressure at 300K','Vapor Pressure at 300K','Pressure at 310K','Vapor Pressure at 310K','Fontsize',14)
grid on
