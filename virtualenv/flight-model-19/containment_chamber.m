%% Establishing Constants
Cv = .50499638; %/coefficient of flow for A = .02 /%
D = 3/16; %/ diameter of outlet, in
D2 = D * .0254; %/diameter of outlet, converted to m
area = (D2/2).^2 * pi; %/ calculates area of outlet
T = 200; %/ temperature of surrounding air
gammaW = 1.33; %/ ratio of specific heats
femtaMassFlow = 80 * 10.^-9; %/mass flow kg/s
Rw = .4615*1000; %/ gas constant for water
V = .05; %/ volume of femta chamber %/ pressure of femta chamber
initialPressure = 0; %/ initial pressure of femta chamber 
tspan = 0:1:150; %Timespan
piraniResolution = 1 * 10.^-3; %Minimum resolution of mks905 micropirani, Pa
P=[];

%% Setting up diff. eq
A= femtaMassFlow * Rw * T / V
B= Rw * T / V * ((Cv * area / sqrt(T)) * sqrt(gammaW/Rw) * ((gammaW + 1) / 2).^(-(gammaW+1)/(2*(gammaW-1))))
C= -A/B

%% Plot pressure vs time
for i = 1:length(tspan)
    P(i) = A/B + C*exp(-B*tspan(i));
end

plot(P);
xlabel('Time (s)');
ylabel('Pressure (Pa)');
title('Pressure vs. Time Inside FEMTA Container');
grid on
grid minor

%/ Check if the sensor can sense the pressure difference from second to second
for i=2:length(P)
    if abs( P(i) - P(i-1)) >= piraniResolution
        fprintf('Success! Pressure differential is %f.6\n',abs(P(i)-P(i-1)));
        fprintf('Timestamp: %f\n',i);
    end
end

% for i = 1:length(P)
%     fprintf('%f\n',P(i));
% 
% end

