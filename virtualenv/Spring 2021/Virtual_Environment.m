clear
clc
close all
%https://multimedia.3m.com/mws/media/199819O/3m-novec-7200-engineered-fluid-en.pdf

%{
Created on Tue Feb 09 2020
    Virtual Environment Simulation of FEMTA propellant management system 
    for testing of flight computer and verification of experimental
    process. Simulation accepts flight conditions and returns relevant data
    for sensor testing under six conditions.
    
    1.) From take-off to the point of experiment startup. This data is
    largely assumed ambient and consists of pressure and temperature. No
    experimental changes contribute to these values. 
    Variable condition for this phase: flightCond=1
    2.) From experiment start to end. There are four possible valve states
    that make up this condition.
    Variable condition for this phase: flightCond=2
        a.) Variable condition for both flow solenoids open and vent 
            solenoid closed (normal): simCond = 1
        b.) Variable condition for one flow solenoid open and vent 
            solenoid closed: simCond = 2
        c.) Variable condition for both fow solenoids open and vent 
            solenoid open: simCond = 3
        d.) Variable condition for one flow solenoid open and vent 
            solenoid open: simCond = 4
    3.) From end of experiment through descent to landing. This data is
    largely assumed ambient and consists of pressure and temperature. No
    experimental changes contribute to these values. 

@author: Alan J
@co-author: Daniel Q
%}

%"tank" suffix denotes propellant tank values
%"CC" suffix denotes collection chamber values
%"O" suffix denotes orifice values
%"pipe" suffix denotes piping network values

%CONSTANT AND INITIAL VALUES
global R kB N_a;
R=8.3145;                 %Universal Gas Constant [kJ/kmol-K]
kB=1.380649e-23;          %Boltzmann constant [J/K]
N_a=6.02214e23;           %Avogadros number [# of particles/mol]

%Air Properties
gammaAir=1.4;             %ratio of specific heats for air
MW_Air=29;                %Molecular Weight of Air [kg/kmol]
R_air=R/MW_Air;           %Specific Gas Constant for Air [kJ/kg-K]
Cp_Air=1.005;             %Specific Heat of Air [kJ/kg-K]
 
%Water Properties
rho_water=997;            %density of liquid water [kg/m^3]
m_H2O=2.988e-26;          %mass of one water molecule [kg]
MW_Water=18.0135;         %Molecular weight of water [kg/kmol]
Cp_water_liquid=4;        %specific heat of liquid water [kJ/kg-K]
Cp_water_vapor=2;         %specific heat of water vapor [kJ/kg-K]
Ce_water=1;               %Evaporation Coefficient of water
Cc_water=0.75;            %Condensation Coefficient of water

%HFE Properties
MW_HFE=250;               %Molecular Weight of HFE [kg/kmol]
R_HFE=R/MW_HFE;           %Specific Gas Constant for HFE [kJ/kg-K]
h_evap_HFE=0.0308/MW_HFE; %Heat of Vaporization of HFE [kJ/kg]
Cp_HFEliquid=1.172303;    %Specific Heat of HFE liquid [kJ/kg-K]
m_HFE=MW_HFE/(N_a*1000);  %mass of one HFE molecule [kg]
Ce_HFE=1;                 %Evaporation Coefficient of HFE
Cc_HFE=1;                 %Condensation Coefficient of HFE

%Propellant Tank
P0_tank=101325;                                           %Initial pressure in prop tank [Pa] (1 atm)
T0_tank=300;                                              %Initial temperature in prop tank [K]
V_tank=(14e-6)*2;                                         %Total volume of both prop tanks [m^3] (30mL) (CONSTANT)
volWater0_tank=9.5e-6;                                    %initial volume of water in one prop tank [m^3] (14mL)
volHFE_liquid0_tank=2.294e-6;                             %initial volume of HFE in both prop tanks [m^3] (10mL)
volAir0=V_tank-(2*volWater0_tank)-volHFE_liquid0_tank;    %Initial volume of air in both prop tanks [m^3]
n_Air=(P0_tank*volAir0)/(R*T0_tank);                      %Number of moles of air in both prop tanks (CONSTANT)
A_HFE=3.167e-5;                                           %Area from which HFE condenses and evaporates [m^2]

%Collection Chamber
P0_CC=0;                  %initial pressure in collection chamber [Pa]
T0_CC=300;                %initial temperature in collection chamber [K]
V_CC=542.248e-6;          %volume of CC [m^3] (227.75mL)
ventSolenoidDiam=0.1e-3;  %Daimeter of vent solenoid [m]

%Piping network
D_pipe=(1/8)/39.37;       %Pipe diameter [m] (CONSTANT)
A=pi*((D_pipe/2)^2);      %Area of pipe cross section [m^2] (CONSTANT)

%Orifice
D_O=0.000127;             %Orifice diameter [m]
A_O=pi*((D_O/2).^2);      %Area of orifice [m^2]
Beta=D_O/D_pipe;          %Ratio of orifice to pipe diameter
CD_orifice=0.6;           %discharge coefficient of orifice

%VARIABLES
volWater_tank=volWater0_tank;                      %Initial volume of water in one Prop Tank [m^3]
volWater_CC=0;                                     %Initial volume of water in Collection Chamber [m^3]
volHFE_liquid=volHFE_liquid0_tank;                 %Volume of HFE Liquid in Prop Tank [m^3]
tankPress=P0_tank;                                 %Pressure in Prop Tank [Pa]
CCPress=P0_CC;                                     %Pressure in Collection Chamber [Pa]
tankTempGas=T0_tank;                               %Temperature of gas in Prop Tank [K]
tankTempLiquid=T0_tank;                            %Temperature of HFE Liquid in Prop Tank [K]
CCTempGas=T0_CC;                                   %Temperature of gas in CC [K]
CCTempLiquid=T0_CC;                                %Temperature of Water Liquid in CC [K]
m_HFE_vapor=0;                                     %Initial mass of HFE vapor in prop tank
m_HFE_liquid=volHFE_liquid0_tank*nvcRho(T0_tank);  %Initial mass of HFE liquid in prop tank
m_water_vapor=0;                                   %Initial mass of water vapor in CC
m_water_liquid=0;                                  %Initial mass of water liquid in CC
n_Gas=n_Air;                                       %Initial moles of gas in tank
volGas=V_tank-(2*volWater_tank)-volHFE_liquid;     %Initial volume of gas in tank [m^3]
volWater_shut=0;                                   %Initial volume of water in tank if one flow solenoid does not open [m^3]
time=0;
count=1;
dt=1e-3;    %timestep [s]

%ARRAY INITIALIZATION
tankVolWater_array=[];
tankPress_array=[];
tankTempGas_array=[];
tankTempLiquid_array=[];
CCVolWater_array=[];
CCPress_array=[];
CCTempGas_array=[];
CCTempLiquid_array=[];
PvapHFE_array=[];
PvapWater_array=[];
QHFE_array=[];
Qwater_array=[];
time_array=[];
m_HFE_transfer_array=[];
m_HFE_vapor_array=[];
m_HFE_liquid_array=[];
tankVolGas_array=[];
n_Gas_array=[];
m_water_vapor_array=[];
m_water_liquid_array=[];
m_water_transfer_array=[];
m_water_unaltered_array=[];
h_evap_water_array=[];
m_HFE_unaltered_array=[];
volGas_array=[];
nGas_array=[];
m_HFE_total_array=[];
m_water_total_array=[];

%COMPUTER SYSTEM STATUS INPUT
%The Sim Box will generate the conditions under which the experiment will
%operate. There are six possible conditions.
%
%Pre-Experiment; ascending to zero-g: flightCond = 1
%During Experiment: flightCond = 2
%   Both flow solenoids open and vent solenoid closed (normal): simCond = 1
%   One flow solenoid open and vent solenoid closed: simCond = 2
%   Both fow solenoids open and vent solenoid open: simCond = 3
%   One flow solenoid open and vent solenoid open: simCond = 4
%Post-Experiment; descending to landing: flightCond = 3
flightCond=2;
simCond=3;

%EXECUTABLE LOOP
if flightCond==1
elseif flightCond==2
        while volWater_tank>0
        
        %Condition if solenoid valve shuts mid experiment
        if (simCond==2 || simCond==4) && (volWater_shut < volWater_tank)
            volWater_shut=volWater_tank;
        end

        %Volumetric Flow Rate of Liquid Water Propellant through one orifice[m^3/s] 
        flo_water=CD_orifice*A_O*sqrt(2*(tankPress-CCPress)/(rho_water*(1-(Beta^4))));

        %Volume lost through vent orifice if open [m^3] (assuming only water
        %vapor lost through vent solenoid)
        if m_water_vapor~=0
            wvRho=m_water_vapor/(V_CC-volWater_CC);
            m_lost=mDotThruOrifice(CCPress,0.000001,wvRho,1.33,1,ventSolenoidDiam)*dt;
            vol_lost=m_lost/wvRho;
        else
            vol_lost=0;
        end

        %Update of Liquid Water Propellant Volumes in Tank & CC [m^3]
        volWater_tank=volWater_tank-(flo_water*dt);
        if simCond==1
            volWater_CC=volWater_CC+(2*flo_water*dt);
        elseif simCond==2
            volWater_CC=volWater_CC+(flo_water*dt);
        elseif simCond==3
            volWater_CC=volWater_CC+(2*flo_water*dt)-vol_lost;
        elseif simCond==4
            volWater_CC=volWater_CC+(flo_water*dt)-vol_lost;
        else
            fprintf("\nError: Invalid Sim Condition\n\n")
            break
        end


        %PROPELLANT TANK
            %Vapor Pressure of HFE [Pa]
            Pvap_HFE=nvcVP(tankTempLiquid);

            %Density of HFE liquid [kg/m^3]
            rho_HFE=nvcRho(tankTempLiquid);

            %Temperature and Pressure of Gas [K, Pa]
            tankPress=(n_Gas*R*tankTempGas)/volGas;
            %tankTempGas=T0_tank*((volAir0/volGas)^(gammaAir-1));

            %Amount of HFE either condensing or evaporating
            m_HFE_transfer=HerKnu(Pvap_HFE,tankTempLiquid,tankTempGas,tankPress,m_HFE,A_HFE,Ce_HFE,Cc_HFE)*dt;
            m_HFE_unaltered=m_HFE_transfer;

            %No mass evaporating if no liquid and no mass condensing if no vapor
            if (tankPress>Pvap_HFE && m_HFE_transfer>0) || (m_HFE_vapor==0 && m_HFE_transfer<0) || (m_HFE_liquid==0 && m_HFE_transfer>0)
                m_HFE_transfer=0;
            end

            %Update amount of HFE liquid and vapor
            m_HFE_liquid=m_HFE_liquid-m_HFE_transfer;
            m_HFE_vapor=m_HFE_transfer+m_HFE_vapor;

            %Mass of vapor and Liquid is zero if negative
            if m_HFE_liquid<0
                m_HFE_liquid=0;
            elseif m_HFE_vapor<0
                m_HFE_vapor=0;
            end

            %Update moles/volume of HFE vapor/liquid 
            vol_HFE_liquid=m_HFE_liquid/rho_HFE;
            n_HFE_vapor=m_HFE_vapor/MW_HFE;

            %Update total amount of Gas (Air + HFE)
            n_Gas=n_HFE_vapor+n_Air;
            if simCond==1 || simCond==3
                volGas=V_tank-(2*volWater_tank)-vol_HFE_liquid;
            else
                volGas=V_tank-(volWater_shut+volWater_tank)-vol_HFE_liquid;
            end

            %Temperatrue Update [K]
            Q_HFE=m_HFE_transfer*h_evap_HFE;

            if m_HFE_liquid>0
                tankTempLiquid=(-Q_HFE/(m_HFE_liquid*Cp_HFEliquid))+tankTempLiquid;
            end
            
            %Check total mass of HFE in Prop Tank [kg]
            m_HFE_total=m_HFE_liquid+m_HFE_vapor;

        %COLLECTION CHAMBER
            %Surface Area of Collected Water [m^2]
            r=((3*volWater_CC)/(4*pi))^(1/3);
            A_water=4*pi*(r^2);

            %Vapor Pressure of Water [Pa]
            Pvap_water=waterVP(CCTempLiquid);

            %Evaporation Heat of Water [kJ/kg]
            h_evap_water=waterHV(CCTempLiquid)/1000;

            %Mass of water either evaporating or condensing at current timestep [kg]
            m_water_transfer=HerKnu(Pvap_water,CCTempLiquid,CCTempGas,CCPress,m_H2O,A_water,Ce_water,Cc_water)*dt;
            m_water_unaltered=m_water_transfer;

            %No mass evaporating if no liquid and no mass condensing if no vapor
            if (CCPress>Pvap_water && m_water_transfer>0) || (m_water_vapor==0 && m_water_transfer<0) || (m_water_liquid==0 && m_water_transfer>0)
                m_water_transfer=0;
            end

            %Total mass of water vapor and liquid at current time [kg]
            m_water_vapor=m_water_transfer+m_water_vapor;
            m_water_liquid=(volWater_CC*rho_water)-m_water_vapor;

            %Mass of Vapor and Liquid is zero if negative
            if m_water_liquid<0
                m_water_liquid=0;
            elseif m_water_vapor<0
                m_water_vapor=0;
            end

            %Total number of moles of water vapor
            n_water_vapor=m_water_vapor/MW_Water;

            %Pressure Update [Pa]
            CCPress=(n_water_vapor*R*CCTempGas)/(V_CC-volWater_CC);

            %Temperature Update [K]
            if CCTempLiquid>273
                Q_water=m_water_transfer*h_evap_water;
            else
                Q_water=0;
            end

            if m_water_liquid>0
                CCTempLiquid=(-Q_water/(m_water_liquid*Cp_water_liquid))+CCTempLiquid;
            end

        %Update total mass of water in the system [kg]
        if simCond==1 || simCond==3
            m_water_total=m_water_vapor+m_water_liquid+(2*volWater_tank*rho_water);
        else
            m_water_total=m_water_vapor+m_water_liquid+((volWater_shut+volWater_tank)*rho_water);
        end

        %ARRAY UPDATE
            %Propellant Tank
            tankVolWater_array(count)=volWater_tank;
            tankPress_array(count)=tankPress;
            tankTempGas_array(count)=tankTempGas;
            tankTempLiquid_array(count)=tankTempLiquid;
            tankVolGas_array(count)=volGas;
            m_HFE_vapor_array(count)=m_HFE_vapor;
            m_HFE_liquid_array(count)=m_HFE_liquid;
            n_Gas_array(count)=n_Gas;
            m_HFE_unaltered_array(count)=m_HFE_unaltered;
            volGas_array(count)=volGas;
            nGas_array(count)=n_Gas;
            m_HFE_total_array(count)=m_HFE_total;

            %Collection Chamber
            CCVolWater_array(count)=volWater_CC;
            CCPress_array(count)=CCPress;
            CCTempGas_array(count)=CCTempGas;
            CCTempLiquid_array(count)=CCTempLiquid;
            m_water_vapor_array(count)=m_water_vapor;
            m_water_liquid_array(count)=m_water_liquid;
            m_water_transfer_array(count)=m_water_transfer;
            m_water_unaltered_array(count)=m_water_unaltered;
            h_evap_water_array(count)=h_evap_water;
            m_water_total_array(count)=m_water_total;

            %Miscellaneous
            PvapHFE_array(count)=Pvap_HFE;
            QHFE_array(count)=Q_HFE;
            PvapWater_array(count)=Pvap_water;
            Qwater_array(count)=flo_water;
            time_array(count)=time;
            m_HFE_transfer_array(count)=m_HFE_transfer;

            time=time+dt;
            count=count+1; 

            %Test Condition
            if time>75.4030
                simCond=1;
            end
        end
elseif flightCond==3
else
    fprintf("\nError: Invalid Flight Condition\n\n")
end

%Data manipulation
time_array=time_array';
m_water_total_array=m_water_total_array';
m_HFE_total_array=m_HFE_total_array';
volGas_array=volGas_array';
nGas_array=nGas_array';
tankVolWater_array=tankVolWater_array';
tankPress_array=tankPress_array';
tankTempGas_array=tankTempGas_array';
tankTempLiquid_array=tankTempLiquid_array';
tankVolGas_array=tankVolGas_array';
m_HFE_vapor_array=m_HFE_vapor_array';
m_HFE_liquid_array=m_HFE_liquid_array';
n_Gas_array=n_Gas_array';
m_HFE_unaltered_array=m_HFE_unaltered_array';
CCVolWater_array=CCVolWater_array';
CCPress_array=CCPress_array';
CCTempGas_array=CCTempGas_array';
CCTempLiquid_array=CCTempLiquid_array';
m_water_vapor_array=m_water_vapor_array';
m_water_liquid_array=m_water_liquid_array';
m_water_transfer_array=m_water_transfer_array';
m_water_unaltered_array=m_water_unaltered_array';
h_evap_water_array=h_evap_water_array';
PvapHFE_array=PvapHFE_array';
QHFE_array=QHFE_array';
PvapWater_array=PvapWater_array';
Qwater_array=Qwater_array';
m_HFE_transfer_array=m_HFE_transfer_array';


%Plots
figure(1)
plot(time_array,tankVolWater_array*10^6,'Linewidth',3)
hold on
plot(time_array,CCVolWater_array*10^6,'Linewidth',3)
xlabel("Time [s]",'Fontsize',17)
ylabel("Volume of Water in Prop Tank and CC [mL]",'Fontsize',17)
legend("One Prop Tank","CC",'Fontsize',15)
title("Volumes of Collection Chamber and Propellant Tank",'Fontsize',22)

figure(2)
plot(time_array,tankTempGas_array)
hold on
plot(time_array,CCTempGas_array)
plot(time_array,CCTempLiquid_array)
plot(time_array,tankTempLiquid_array)
xlabel("time [s]",'Fontsize',17)
ylabel("Temperature [K]",'Fontsize',17)
legend("Temp of gas in tank","Temp of gas in CC","temp of liquid in CC","Temp of liquid in tank",'Fontsize',17)
title("Temperature values for Duration of Experiment",'Fontsize',22)

figure(3)
plot(time_array,tankPress_array/1000,'Linewidth',3)
hold on
plot(time_array,PvapHFE_array/1000,'Linewidth',3)
xlabel("Time [s]",'Fontsize',17)
ylabel("Pressure of Prop Tank [kPa]",'Fontsize',17)
legend("Tank Pressure","Vapor Pressure of HFE",'Fontsize',17)
title("Pressures in Propellant Tank",'Fontsize',22)

figure(4)
plot(time_array,CCPress_array/1000,'Linewidth',3)
hold on
plot(time_array,PvapWater_array/1000,'Linewidth',3)
xlabel("Time [s]",'Fontsize',17)
ylabel("Pressure of Collection Chamber [kPa]",'Fontsize',17)
legend("CC Pressure","Vapor Pressure of Water",'Fontsize',17)
title("Pressures in Collection Chamber",'Fontsize',22)

figure(5)
plot(time_array,m_HFE_vapor_array*1000)
hold on
plot(time_array,m_HFE_liquid_array*1000)
plot(time_array,m_HFE_total_array*1000)
xlabel("Time [s]")
ylabel("Mass of HFE [g]")
title("Amount of HFE in Each State")
legend("HFE Vapor","HFE Liquid","HFE Total")

figure(6)
plot(time_array,m_water_vapor_array*1000)
hold on
plot(time_array,m_water_liquid_array*1000)
plot(time_array,m_water_total_array*1000)
plot(time_array,tankVolWater_array*rho_water*1000)
xlabel("Time [s]")
ylabel("Mass of Water [g]")
title("Amount of Water in Each State")
legend("Water Vapor","Water Liquid","Water Total","Water in One Prop Tank")

% figure(7)
% plot(time_array,m_water_unaltered_array)

%Vapor Pressure of NV 7100 (T in K, vp in Pa)
function vp = nvcVP(T)
    vp = exp(22.415 - 3641.9 * (1/T));
end

%Density of NV 7100 Liquid (T in K)
function dD = nvcRho(T)
    dD = 1.5383 - 0.002269*(T-273.15);
    dD = dD / (1000 * 0.000001); %kg/m^3
end

%Vapor Pressure of Water (T in K, vp in Pa)
function vp = waterVP(T)
    vp=10.^(8.07131-(1730.63./(233.426+(T-273))));        %pressure in mmHg
    vp=vp*133;                                            %conversion to Pa from mmHg
end

%Mass transfer from gas to liquid [kg/s] (HERTZ-KNUDSEN EQUATION)
%Negative denotes vapor to liquid (condensation), Positive denotes liquid to vapor
%(evaporation)
function m_transfer=HerKnu(Ps,T_liquid,T_vapor,Pg,m,A,C_evap,C_cond)
    global kB;
    m_transfer=A*(sqrt(m/(2*pi*kB))*((C_evap*(Ps/sqrt(T_liquid)))-(C_cond*(Pg/sqrt(T_vapor)))));
end

%Heat of Vaporization of Water [J/kg]
function hv=waterHV(T)
    Hvs=[2500.9 2496.2 2491.4 2477.2 2467.7 2458.3 2453.5 2441.7 2429.8 2420.3 2406 2396.4 2381.9 2372.3 2357.7 2333 2308 2282.5 2266.9 2256.4 2229.6 2202.1 2144.3 2082 2014.2 1939.7 1857.4 1765.4 1661.6 1543 1404.6 1238.4 1027.3 719.8];
    Ts=[0.00 2 4 10 14 18 20 25 30 34 40 44 50 54 60 70 80 90 96 100 110 120 140 160 180 200 220 240 260 280 300 320 340 360]+273;
    hv=interp1(Ts, Hvs, T, 'linear').*1000;
end

%mDotThruOrifice calculates mass flow (mDot) in kg/s through an orifice
%given pressures and working fluid properties
function [mDot] = mDotThruOrifice(in1,in2,in3,in4,in5,in6)
%   Variable description:
%       in1: pressure upstream (Pa)
%       in2: pressure downstream (Pa)
%       in3: fluid density (kg/m^3)
%       in4: gamma = (fluidCP/fluidCV)
%       in5: orifice discharge coefficient
%       in6: orifice diameter (m)
%
%   Note on in1 and in2 variables:
%       The sign convention used in this function assumes normal (positive)
%       fluid movement from the in1 region to the in2 region (in1 is
%       upstream, in2 is downstream). However, these variables may be 
%       reversed so that there is flow from in2 to in1, however, the
%       resulting mDot will be negative.
%   
    if in1 < in2
        downP = in1;
        upP = in2;
        directionSign = -1;
    else
        upP = in1;
        downP = in2;
        directionSign = 1;
    end
    rho = in3;
    gamma = in4;
    outletCD = in5;
    outletDia = in6;
    outletArea = pi*(outletDia/2)^2;
    %fprintf("%0.3f, %0.3f, %0.6f, ", upP, downP, rho);
    %
    criticalP = upP * (2/(gamma+1))^(gamma/(gamma-1));
    %fprintf("%0.3f, %0.3f, ", criticalP, downP);
    if(downP < criticalP)
        %Choked
        %fprintf("C\n");
        r = downP/criticalP;
        r = (2/(gamma+1))^(gamma/(gamma-1));
        mDot = outletCD*outletArea*sqrt(upP*rho*(2*gamma/(gamma-1))*r^(2/gamma)*(1-r^((gamma-1)/gamma))); %kg/s
        %fprintf("%0.20f, ", rho);
    else
        %Subsonic
        %fprintf("S\n");
        r = downP/upP;
        mDot = outletCD*outletArea*sqrt(upP*rho*(2*gamma/(gamma-1))*r^(2/gamma)*(1-r^((gamma-1)/gamma))); %kg/s
    end
    %fprintf("\n");
    mDot = mDot * directionSign; %Corrects sign on mDot to follow stated convention above
end

