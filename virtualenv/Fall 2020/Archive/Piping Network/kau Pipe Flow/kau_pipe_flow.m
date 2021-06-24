clear
clc
%VARIABLES
L=0.05;                 %Length of horizontal pipe [m]
mu=0.0006527;           %Dynamic Viscosity [Pa*s] (@ 40C)
nu=1.787e-6;            %Kinematic Viscosity of water [m^2/s] (@ 0C)
V_c=10;                 %Centerline velocity [m/s]
P1=1;                   %Pressure at propellant tank exit [Pa]
R=(1/8)*0.5*(1/39.37);  %Pipe Radius [m]
D=2*R;                  %Pipe Diameter [m]
A=pi*(R^2);             %Pipe Area [m^2]
V0=14;                  %Initial volume of propellant tank [m^3]
Vf=0;                   %Final volume of propellant tank [m^3]
tspan=linspace(0,150);  %time that water flows out of tank
rho=997;                %density of water [kg/m^3]
epsilon=0.01e-3;        %average pipe roughness [m] (for plastic)
K=0.3;                  %Minor Loss Factor due to 90 degree bend
g=9.8;                  %Acceleration of gravity [m/s^2]

%Volumetric Flow Rate
V_t=(((Vf-V0)/tspan(length(tspan)))*tspan)+V0;
dvdt=(-(Vf-V0)/tspan(length(tspan)))/(1e6);

%SECTION 1
del_P1=(1.07*rho*(dvdt^2)*L/(D^5))*(log((epsilon/D/3.7)+(4.62*((nu*D/dvdt)^0.9)))^-2);

%SECTION 2
del_P2=K*rho*((dvdt/A)^2)/2;

%TOTAL
total_del_P=sum([del_P1 del_P2])
