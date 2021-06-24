%Code created by Alan Johnson
clear
clc
%Global Definitions
global m kB rhol N_a MW D A R C_evap C_cond V_CC h_evap h_c;        
m=2.988e-26;        %mass of one water molecule [kg]       
kB=1.380649e-23;    %Boltzmann constant [J/K]      
rhol=997;           %density of liquid water [kg/m^3]    
N_a=6.02214e23;     %Avogadros number [# of particles/mol]
MW=0.0180153;       %Molecular weight of water [kg/mol]
D=(1/8)/39.37;      %Pipe diameter [m]
A=pi*((D/2)^2);     %Area of pipe cross section [m^2]
R=8.314;            %Universal Gas Constant of Water [J/mol-K]
C_evap=0.01;           %Evaporation Coefficient
C_cond=0.01;           %Condensation Coefficient
V_CC=38e-6;         %volume of CC [m^3]
h_evap=2500;        %heat of vaporazitaion of water [kJ/kg]
h_c=4000;           %convective heat transfer coefficient of water [kW/m^2-K]

%Iteration Starters
dt=0.003;           %timestep [s]
Pg=0.00001;         %initial Pressure in CC
dT=274;             %initial temp [K]
m_liquid_last=0;    %initial mass of liquid water in CC
m_vapor_last=0;     %initial mass of water vapor in CC
fluid_density=0;    %initial density of fluid in CC
totalTime=150;      %duration of experiment

%Storage Vectors
N=round(totalTime/dt)+1;
count=1;
CCPressure=zeros(N,1);
massLiquid_array=zeros(N,1);
massVapor_array=zeros(N,1);
volumeVapor_array=zeros(N,1);
Pg_array=zeros(N,1);
HerKnu_array=zeros(N,1);
orifice_array=zeros(N,1);
temp_array=zeros(N,1);

%ANTOINE EQUATION CONSTANTS (valid for 274-373 K) (generates P in mmHg)
Z=8.07131;
B=1730.63;
C=233.426;
%ANTOINE EQUATION CONSTANTS (valid for 374-647 K) (generates P in mmHg)
E=8.14019;
F=1810.94;
G=244.485;

for simTime=0:dt:totalTime
    
    %Vapor Pressure of Water
    if dT<=373
        Pv=10.^(Z-(B./(C+(dT-273))));        %interface pressure in mmHg
        Pv=Pv*133;                           %conversion to Pa from mmHg
    elseif dT>374
        Pv=10.^(E-(F./(G+(dT-273))));        %interface pressure in mmHg
        Pv=Pv*133;                           %conversion to Pa from mmHg
    end
    
    %mass of liquid water in CC [kg]
    mdot_liquid=((14/totalTime)/(1e6))*rhol;     %Assuming constant linear mass flow
    m_liquid=(mdot_liquid*dt)+m_liquid_last;
    massLiquid_array(count)=m_liquid;
    m_liquid_last=m_liquid;
    
    %volume of liquid in CC [m^3]
    V_liquid=m_liquid/rhol;
    
    %mass of water vapor in CC [kg]
    m_vapor=((HerKnu(Pv,dT,Pg)-mDotThruOrifice(Pg,0,fluid_density,1.33,0.01,0.0001))*dt)+m_vapor_last;
    massVapor_array(count)=m_vapor;
    m_vapor_last=m_vapor;
    HerKnu_array(count)=HerKnu(Pv,dT,Pg);
    orifice_array(count)=mDotThruOrifice(Pg,0,fluid_density,1.33,0.01,0.01);
    
    %volume of gas in CC [m^3]
    V_vapor=V_CC-V_liquid;
    volumeVapor_array(count)=V_vapor;
    fluid_density=m_vapor/V_vapor;
    
    %Pressure in CC at dt [Pa]
    Pg=(m_vapor*R*dT)/(MW*V_vapor);
    
    %Temperature Change due to Evaportion
    temp_array(count)=dT;
    dT_last=dT;
    dT=((h_evap*HerKnu(Pv,dT,Pg))/(h_c*A))+dT_last;
    
    %Storage Vector Fill
    CCPressure(count)=Pg;
    count=count+1; 
end

plot((0:dt:totalTime),CCPressure,'Linewidth',3)
title("Collection Chamber Pressure for Duration of Experiment",'Fontsize',22)
ylabel("Pressure [Pa]",'Fontsize',17)
xlabel("Time [s]",'Fontsize',17)
grid on
%set(gca,'Color','k')

[massLiquid_array massVapor_array];
HerKnu_array;
temp_array;
CCPressure(1:5);
Pv;
%Mass transfer from gas to liquid [kg/s]
%Negative denotes vapor to liquid (condensation), Positive denotes liquid to vapor
%(evaporation)
function m_transfer=HerKnu(Pv,dT,Pg)
    global m kB A C_evap C_cond;
    m_transfer=A*(sqrt(m/(2*pi*kB))*((C_evap*(Pv/sqrt(dT)))-(C_cond*(Pg/sqrt(dT)))));
end

%Radius of assumed sphere of collected water [m]
% function dr=radius(dmdt)
%     dr=dr+0.001;
% end

%Mass collection rate of water on surface of sphere in collection chamber
%[kg/s]
% function dmdt=mass_collect(dr,dP,dT)
%     A=Area(dr);
%     dmdt=(A*MW*alpha*dP)/(N_a*sqrt(2*pi*m*kB*dT));
% end

%Temperature at given time in collection chamber
% function dT=temp(dt);
%     dT=
% end

%Surface Area of Collection Geometry
% function A=Area(dr)
%     A=4*pi*(dr^2);
% end

%Flow rate of water into collection chamber by mass
% function mdot=mass_flo(P_diff)
%     mdot=rho*(pi*((D/2)^2))*sqrt((2*P_diff)/rho);
% end
