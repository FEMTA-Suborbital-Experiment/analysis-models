%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This algorithm determines the pressure difference along a pipe of given
%length and under the given conditions.
%This method is for the assumption of a horizontal, level pipe with no
%obstructions to the flow. The length and dynamic viscosity are currently
%arbitrary. Volumetric flow rate is assumed to be constant and linear for
%the duration of the 430s time span. Laminar flow is assumed with a parablic
%velocity profile. See the page "Pipe Flow Equations" for derivation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%VARIABLES
L=0.05;                 %Length of horizontal pipe [m]
mu=8.9e-4;              %Dynamic Viscosity [Pa*s] (@ 40C)
V_c=10;                 %Centerline velocity [m/s]
P1=1;                   %Pressure at propellant tank exit [Pa]
R=(1/8)*0.5*(1/39.37);  %Pipe Radius [m]
V0=0.014;               %Initial volume of propellant tank [m^3]
Vf=0;                   %Final volume of propellant tank [m^3]
tspan=linspace(0,10);   %time that water flows out of tank

% %USING FLOW VELOCITY
% del_P1=(4*L*mu*V_c)/(R^2);
% P2=del_P1+P1;

%USING VOLUMETRIC FLOW RATE
V_t=(((Vf-V0)/tspan(length(tspan)))*tspan)+V0;
dvdt=((Vf-V0)/tspan(length(tspan)))/(1e6);
del_P=double((-dvdt*128*mu*L)/(pi*((2*R)^4)))
speed=-dvdt/(pi*(R^2));
