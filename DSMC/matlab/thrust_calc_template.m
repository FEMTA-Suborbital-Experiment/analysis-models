clear
clc
%Load data from SPARTA
femtaNrho=load('femtaNrho_timestep800.out');
femtaPressTemp=load('femtaPressTemp_timestep800.out');

%Define Constants
Na=6.0221409e23;      %avogadros number [# of particles/mol]
MW=0.0180153;         %molar weight of water [kg/mol]
fnum=1e5;
d=1.25/2.54;          %cathode diameter at mouth [cm]                         

%Grid cell center locations
xc=femtaNrho(:,1);
yc=femtaNrho(:,2);
zc=femtaNrho(:,3);


%Number Density, Pressure, and Velcoity for entire sim box
total_count=femtaNrho(:,4);
total_nrho=femtaNrho(:,8);
total_pressure=femtaPressTemp(:,4);
total_xvel=femtaNrho(:,5);
total_yvel=femtaNrho(:,6);
total_zvel=femtaNrho(:,7);

%define vectors to fill
flux_nrho=zeros(1,length(zc));
flux_pressure=zeros(1,length(zc));
flux_velocity=zeros(1,length(zc));
count=1;

%Defines number density, pressure, and velocity for all x and y at z
%nearest to femta exit
for i=1:length(zc)
    if zc(i)==-0.00113419
        flux_nrho(count)=total_nrho(i);
        flux_pressure(count)=total_pressure(i);
        flux_velocity(count)=sqrt((total_xvel(i)^2)+(total_yvel(i)^2)+(total_zvel(i)^2));
        count=count+1;
    end
end

%Area of each individual grid cell
A=((0.021*0.031)*5)/count;      %area of each cell in the flux region[m^2]

%Convert number density to density
flux_rho=((MW/Na)*flux_nrho);

%Strip trailing zeros from vectors
rho2=find(flux_rho,1,'last');
flux_rho=flux_rho(1:rho2);
press2=find(flux_pressure,1,'last');
flux_pressure=flux_pressure(1:press2);
vel2=find(flux_velocity,1,'last');
flux_velocity=flux_velocity(1:vel2);

%Calculate Thrust
Thrust=(flux_rho.*(flux_velocity.^2).*A)+(flux_pressure.*A);

%Pressure at cthode mouth
P=(6e-6);     %[Pa]

fprintf('The calculated value of thrust for FEMTA is %.4f microNewtons.\n', sum(Thrust)*(10^6)*4)
fprintf('The average pressure at the cathode mouth is %e Torr-cm.\n', (P/133)*d)

% scatter3(xc,yc,zc)
% xlabel("x")
% ylabel("y")
% zlabel("z")



