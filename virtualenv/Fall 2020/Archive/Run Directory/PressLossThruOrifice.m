%Source
%https://neutrium.net/fluid-flow/calculation-of-flow-through-nozzles-and-orifices/

%Constants
DO=0.0001;  %Orifice diameter [m]
AO=pi*((DO/2).^2);           %Area of orifice [m^2]
D1=(1/8)*(1/39.37);          %Pipe Diameter [m]
Cd=0.8;                      %Discharge Coefficient (thin, sharp-edged orifice plate)
Beta=DO/D1;                  %Ratio of orifice to pipe diameter
Y=1;                         %Expansion coefficient (assuming incompressible flow)
rho=997;                     %density of liquid water [kg/m^3]
Vf=0;                        %final volume of liquid in propellant tank [mL]
V0=14;                       %initial volume of liquid in propellant tank [mL]
tspan=0:0.03:199.97;       %time of experiment [s]

%Volumetric Flow rate through orifice
dvdt=volumetric_flow_rate;   %from tank pressure

%Pressure Loss
P1=1;                %Upstream Pressure [Pa]
delP_inOrifice=0.5*rho*(1-(Beta.^4)).*((dvdt./(Cd*AO*Y)).^2);
%fprintf('The pressure loss though the orifice is %.4f Pa.\n', delP)

plot(tspan,delP_inOrifice/1000,'Linewidth',3)
ax=gca;
ax.FontSize=18;
xlabel('Time [s]','Fontsize',20)
ylabel('Pressure Change [kPa]','Fontsize',20)
title('Pressure Change Across an Orifice','Fontsize',22)
grid on


