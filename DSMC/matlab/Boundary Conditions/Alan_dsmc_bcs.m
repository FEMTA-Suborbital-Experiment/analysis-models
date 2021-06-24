clear
clc
%VARIABLE DECLARATION
R=461.52;                   %Specific Gas Constant [J/kg-K]
avg=6.022e23;               %Avogadro's Number
gamma=1.257;                %ratio of specific heats
eta=.866;                   %enthalpy efficiency of nozzle
Ti=[295:0.1:315];           %temperature range [K]
Ae = 0.000000008;           %exit area for FEMTA
MW_H2O = 0.01801524;        %molar mass of water [kg/kmol]
%ANTOINE EQUATION CONSTANTS (valid for 293-343 K) (generates P in bar)
A=6.20963;
B=2354.731;
C=7.559;
%INTERFACE PRESSURE (vapr pressure of water at Ti)
Pi=10.^(A-(B./(C+Ti)));                     %interface pressure in bar
Pi=Pi*100000;                               %conversion to Pa from bar 
%ITERATIVE PROCESS TO CALCULATE Texit (Te in K) (cpw calculates specific heat)
cp1=cpw(Ti);
TeInt1=Ti.*(1-eta.*(R./(2.*cp1-R)));
cp2=cpw(TeInt1); 
TeInt2=Ti.*(1-eta*(R./(2*cp2-R)));          
cp3=cpw(TeInt2);
Te=Ti.*(1-eta*(R./(2*cp3-R)));          %BC1
%EXIT VELOCITY
cp_Ti=cpw(Ti);                              %Specific heat at interface T
cp_Te=cpw(Te);                              %Specific heat at exit T       
Ue=sqrt(2*cp_Te.*(Ti-Te));                  %Exit velocity (m/s) (BC2)
%EXIT DENSITY       
rho_e=((Te./Ti).^(cp_Te/R)).*(Pi./(R*Te));  %mass density at exit [kg/m^3]
n_e=(rho_e/MW_H2O)*avg;                     %number desnity at exit [#/m^3] (BC3)
%EXIT THRUST
Th=(Ae*Pi.*((Te./Ti).^(cp_Te./R)).*(((Ue.^2)./(R*Te))+1))*(10^6);   %Thrust [N] (BC4)

%PLOTS
subplot(2,2,1)
plot(Ti,Th,'.')
title('Thrust vs. Interface Temperature','Fontsize',17);
xlabel('Interface Temperature [K]','Fontsize',15)
ylabel('Thrust [\muN]','Fontsize',15)
subplot(2,2,2)
plot(Ti,Te,'.');
xlabel('Interface Temperature [K]','Fontsize',15);
ylabel('FEMTA Exit Temperature [K]','Fontsize',15);
title('Exit Temperature vs. Interface Temperature','Fontsize',17);
subplot(2,2,3)
plot(Ti,n_e,'.');
xlabel('Interface Temperature [K]','Fontsize',15);
ylabel('Number Density [particles/m^3]','Fontsize',15);
title('Number Density vs. Interface Temperature','Fontsize',17);
subplot(2,2,4)
plot(Ti,Ue,'.');
xlabel('Interface Temperature [K]','Fontsize',15);
ylabel('Velocity [m/s]','Fontsize',15);
title('Exit Velocity vs. Interface Temperature','Fontsize',17);

fprintf('The average calculated value of thrust for FEMTA is %.4f microN.\n', mean(Th))
fprintf('The average calculated value of exit temperature for FEMTA is %.4f K.\n', mean(Te))
fprintf('The average calculated value of exit velocity for FEMTA is %.4f m/s.\n', mean(Ue))
fprintf('The average calculated value of number density for FEMTA is %e particles/m^3.\n', mean(n_e))




