%% NOTE: This code is still in the works. Specifically, the high dependency
%        of fluid flow rate on the diameter of the orifice is not modeled
%        here

%% Setup
collectionP = 3000; %pa
pipeID = 0.050; %in (inside pipe diameter)
valveCv = 0.032; %imperial units
orificeCv = 0.0012; %imperial units
orificeCd = 0.8; 
orificeD = 0.0071; %in

tankV = 316.254; %cm^3
fluidT = 280.15; %K (working fluid temperature)
rho = 1000; %kg/m^3 (working fluid density)
fluidSG = 1; % (specific gravity working fluid)
hfeT = 280.15; %K (pressurant fluid temperature)
hfeVapP = exp(22.415-3641.9*(1/hfeT)); %Pa (vapor pressure of pressurant at temperature)

mDot = 0.00004; %kg/s (highly dependent on orifice size and pressure drop
fV = mDot/(rho*(pi*(pipeID/2/39.37)^2)); %m/s (flow velocity)
hD = pipeID; %in (hydraulic diameter)

mu = 0.00002414 * 10^(247.8/(fluidT - 140)); %Ns/m^2 (dynamic fluid velocity)
Re = (rho * fV * (hD/39.37)) / mu; %Reynolds number
f = 0; %friction factor
if Re >= 3000 && Re <= 3*10^6 %calculating friction factor from reynolds number
    f = 0.0014 + (0.125/(Re^0.32));
elseif Re < 2100
    f = 16/Re;
else
    fprintf('Reynolds Number out of bounds\n');
end

%% Calculations
%dp across all bends in system
allBendDP = 2*pipe_dp(90,0.5/39.37,0,0.05/39.37,rho,fV,f,f); %90 degree bend
%dp from pipe wall friction
allStraightDP = 7*pipe_dp(0,0,1.0/39.37,0.05/39.37,rho,fV,f,f) + pipe_dp(0,0,2/39.37,0.05/39.37,rho,fV,f,f);

%dp across valve
vFR = (1/rho)*mDot; %m^3/s (volumetric flow rate)
valveDP = 1/(((valveCv/(vFR*15850.323))^2)/fluidSG) * 6894.757; %Pa

%dp across orifice plate (need to size correctly for proper mDot)
orificeDP = 1/(((orificeCv/(vFR*15850.323))^2)/fluidSG) * 6894.757; %Pa;
orificeMDot = (0.8/(sqrt(1-((orificeD/39.37)/(pipeID/39.37))^4)))*1*(pi/4)*((orificeD/39.37)^2)*sqrt(2*rho*orificeDP);

tempOrificeDP = ((4*mDot*sqrt(1-(orificeD/pipeID)^4))/(orificeCd*1*pi*((orificeD/39.37)^2)*sqrt(2*rho)))^2;

flowMeterDP = 0;
totalDP = allBendDP + allStraightDP + valveDP + flowMeterDP + orificeDP;
fprintf('Total dP = %0.3f Pa\n', totalDP);

