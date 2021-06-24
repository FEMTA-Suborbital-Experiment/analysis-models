%{
Created on Fri Sept 22 22:34:15 2020
    Simulated Hardware In The Loop (SHITL) code for the virtual
    environment. The virtual environment's purpose is to test the flight
    computer to ensure proper (and predictable) functionality. This will
    suplement full system testing once all hardware is manufactured and in
    a way serves as a precursor test to work out any issues while hardware
    is being manufactured.
    -> This code take in inputs for valve actuations and outputs all flow
    data throughout the system. Values are most notably pressure and
    temperature, however, other non-measureable values are calculated such
    as masses and partial pressures.
    -> NOVEC Freezing point alert uses the 1 atm freezing point temp.
    -> MUST VERIFY DIAPHRAGM CHARACTERISTICS! (dimensions don't agree!)
    -> Assuming no pressure drop (or negligible). Only driving force for
    water is expanding diaphragm
    -> Modeling surface area transient for liquids?
@author: Justin C (pressure build-up in propellant tank)
@co-author: Alan J (pipe flows (not yet implemented here))
%}

clear
LVTF_ambientP = load('LVTF_ambientP.mat').ans + 0.25;
ambientAlt = load('altitude.mat').h;
flightTime = load('time.mat').t;
livePlot = 0; %bool plot data live (LAGGY!!)
stoData = 1; %bool store data in arrays at every time step

%initial conditions
ctrlChmbrVent = 0;
ctrlTankVent = 0;
ctrlTankRun = 1;
%ambientP = 45000;
[ambientT, ambientP, ambientRho] = StandardAtm(ambientAlt(1));
chmbrTG = 295.6;
chmbrTL = 295.6;
chmbrPNvc = nvcVP(chmbrTL);
chmbrPAir = ambientP-chmbrPNvc;
chmbrML = 0.1*(1e-6)*nvcRho(chmbrTL);
chmbrVL = chmbrML/nvcRho(chmbrTL);
chmbrVG = (5e-7)-chmbrVL;
chmbrSAL = 0.0000149034;
chmbrMGAir = chmbrPAir*chmbrVG/287.05/chmbrTG;
chmbrMGNvc = chmbrPNvc*chmbrVG/33.25/chmbrTG;
chmbrVentDia = 0.001;
tankTG = 295;
tankTL = 295;
tankPH2o = h2oVP(tankTL);
tankPAir = ambientP-tankPH2o;
tankVL = 13.5*(1e-6);
tankVG = 0.5*(1e-6);
tankLRho = 997; %initial density
tankML = tankLRho*tankVL;
tankSAL = 0;
tankMGAir = tankPAir*tankVG/287.05/tankTG;
tankMGH2o = tankPH2o*tankVG/461.52/tankTG;
tankVentDia = 0.001;
tankOutletDia = 0.003175;
tankOutletVFlow = 0; %volumetric flow out of tank outlet (negative: out of tank)
tankOutletMFlow = 0; %mass flow out of tank outlet (negative: out of tank)
tankPBack = backPress(tankPH2o + tankPAir, tankOutletMFlow);

dPLast = 0; %for diaphragm characteristics
KIDiaph = 0; %for diaphragm characteristics

resMG = 10;
resT = 280;
resV = 15*(1e-6);
resPHist = zeros([1,0]);

if livePlot
    clf;
    hold on;
    yyaxis left;
    xlim([0 150]);
    ylim([-5 45]);
    yyaxis right;
    ylim([282 300]);
    grid on;
end
if stoData
    simData.time = zeros([0,0]);
    simData.physics.mass.values = zeros([0,0]);
    simData.physics.mass.labels = [];
    simData.physics.pressure.values = zeros([0,0]);
    simData.physics.pressure.labels = [];
    simData.physics.temperature.values = zeros([0,0]);
    simData.physics.temperature.labels = [];
    simData.physics.volume.values = zeros([0,0]);
    simData.physics.volume.labels = [];
    simData.events.values = zeros([0,2]);
    simData.events.labels = ["Time", "Event"];
end

dt = 0.03;
ffNvcVent = 0.02;
ffH2oVent = 0.005;
ffH2oOutlet = 0.005;
ffDiaph = 0.000000001;
simData.events.values = [simData.events.values; [0, "Simulation start"]];
for simTime = -1.5:dt:149
    %% Control transient
    if(simTime >= 60 && simTime < 60+dt)
        
        ctrlChmbrVent = 1;
        simData.events.values = [simData.events.values; [simTime, sprintf("ctrlChmbrVent State to %d",ctrlChmbrVent)]];
        %{
        ctrlTankVent = 1;
        if(simTime >= 60)
            ctrlChmbrVent = 0;
            if(simTime >= 90)
                ctrlChmbrVent = 1;
            end
        end
        %}
    elseif(simTime >= 80 && simTime < 80+dt)
        ctrlTankVent = 1;
        simData.events.values = [simData.events.values; [simTime, sprintf("ctrlTankVent State to %d",ctrlTankVent)]];
    elseif(simTime >= 150 && simTime < 150+dt)
        ctrlTankRun = 0;
        simData.events.values = [simData.events.values; [simTime, sprintf("ctrlTankRun State to %d",ctrlTankRun)]];
    end
    
    %ambientP = 1000*LVTF_ambientP(find(LVTF_ambientP(:,1) > simTime, 1, 'first'), 2);
    if(simTime < 0)
        [ambientT, ambientP, ambientRho] = StandardAtm(ambientAlt(find(round(flightTime, 4) == round(0.0, 4))));
    else
        [ambientT, ambientP, ambientRho] = StandardAtm(ambientAlt(find(round(flightTime, 4) == round(simTime, 4))));
    end
    %% NOVEC evaporation transient
    nvcEC = 0.001625; %0.01625
    nvcCC = 0.0016314; %0.016314
    if(ctrlChmbrVent == 0)
        nvcEC = 0.0115;
        nvcCC = 0.009;
    end
    %
    mDot = nvcEvap(nvcEC, nvcCC, chmbrSAL, chmbrPNvc, chmbrTG, chmbrTL);
    chmbrML = chmbrML - mDot*dt;
    %fprintf("%0.10f kg\t", chmbrML);
    if(chmbrML < 0)
        mDot = mDot - abs(chmbrML)/dt;
        chmbrML = 0;
    end
    chmbrMGNvc = chmbrMGNvc + mDot*dt;
    chmbrVL = chmbrVL - mDot/nvcRho(chmbrTL)*dt;
    chmbrVG = chmbrVG + mDot/nvcRho(chmbrTL)*dt;
    chmbrPNvc = chmbrMGNvc*33.25*chmbrTG/chmbrVG;
    %% NOVEC thermal transient
    if(chmbrML > 0)
        if(ctrlChmbrVent == 0)
            dT = 0.34 * (112) * -mDot * (1/(1.183*chmbrML));
            dT = dT + 10.0*(0.04*chmbrSAL*(chmbrTL-chmbrTG)*(1/(1.183*chmbrML)));
            chmbrTL = chmbrTL + dT*dt;
        else
            %
            dT = (chmbrML/0.0015)*(112/1000) * -mDot * chmbrSAL * (1/(1.183*chmbrML));
            dT = dT - 0.8*(1*chmbrSAL*((chmbrTL)-chmbrTG)*(1/(1.183*chmbrML)));
            chmbrTL = chmbrTL + dT*dt;
            dT = -0.03*(1*chmbrSAL*(chmbrTG-(chmbrTL))*(1/(1.183*chmbrML)));
            chmbrTG = chmbrTG + dT*dt;
            %
        end
    else
        chmbrML = 0;
        chmbrTL = chmbrTG;
    end
    %fprintf("%0.10f K\n", chmbrTG)
    if(chmbrTL < (-135+273.15))
        simData.events.values = [simData.events.values; [simTime, "NOVEC liquid freezing!"]];
    end
    %% NOVEC vent transient
    if(ctrlChmbrVent == 0)
        mDot = ffNvcVent*mDotThruOrifice(ambientP, chmbrPAir+chmbrPNvc, (chmbrMGAir+chmbrMGNvc)/chmbrVG, (1.4*chmbrMGAir+1.0289*chmbrMGNvc)/(chmbrMGAir+chmbrMGNvc), 0.8, chmbrVentDia);
        chmbrMGAir = chmbrMGAir + mDot*(chmbrMGAir/(chmbrMGAir+chmbrMGNvc))*dt;
        chmbrMGNvc = chmbrMGNvc + mDot*(chmbrMGNvc/(chmbrMGAir+chmbrMGNvc))*dt;
        if(chmbrMGAir < 0)
            chmbrMGAir = 0;
        end
        if(chmbrMGNvc < 0)
            chmbrMGNvc = 0;
        end
        chmbrPAir = chmbrMGAir*287.05*chmbrTG/chmbrVG;
        chmbrPNvc = chmbrMGNvc*33.25*chmbrTG/chmbrVG;
        % Thermal effect
        ventVel = mDot / ((chmbrMGAir+chmbrMGNvc)/chmbrVG*chmbrVentDia);
        dT = 0.5*ventVel^2 * mDot * (1/(1.183*chmbrML));
        dT = dT + (0.08*chmbrSAL*(chmbrTL-chmbrTG)*(1/(1.183*chmbrML)));
        chmbrTG = chmbrTG + dT*dt;
    end
    %% H2O evaporation transient
    %
    if(tankVG > 0)
        if(tankSAL == 0)
            tankSAL = pi*(tankVentDia/2)^2;
        end
        %{
        if(tankMGAir == 0 && tankMGH2o == 0)
            tankPAir = 0;
            tankPH2o = 0;
        end
        %}
        mDot = h2oEvap(0.0015, 0.0015, tankSAL, tankPH2o, tankTG, tankTL);

        tankML = tankML - mDot*dt;
        if(tankML < 0)
            mDot = mDot - abs(tankML)/dt;
            tankML = 0;
        end
        tankVG = tankVG - (tankML/tankLRho-tankVL);
        tankVL = tankVL + (tankML/tankLRho-tankVL);
        if(simTime == 70.11)
            tankVG
        end
        tankLRho = tankML/tankVL;
        tankMGH2o = tankMGH2o + mDot*dt;
        tankPH2o = tankMGH2o*461.52*tankTG/tankVG;
    else
        tankVG = 0;
        tankMGH2o = 0;
        tankPH2o = chmbrPAir+chmbrPNvc;
    end
    %
    %% H2O vent/outlet transient
    if((ctrlTankVent == 0 || ctrlTankRun == 0) && tankVG > 0)
        if(tankSAL == 0)
            tankSAL = pi*(tankVentDia/2)^2;
        end
        if(ctrlTankVent == 0)
            mDot = ffH2oVent*mDotThruOrifice(ambientP, tankPAir+tankPH2o, (tankMGAir+tankMGH2o)/tankVG, (1.4*tankMGAir+1.329*tankMGH2o)/(tankMGAir+tankMGH2o), 0.8, tankVentDia);
        end
        if(ctrlTankRun == 0)
            tankPBack = backPress(tankPAir+tankPH2o, tankOutletMFlow);
            tankOutletMFlow = ffH2oOutlet*mDotThruOrifice(tankPBack, tankPAir+tankPH2o, (tankMGAir+tankMGH2o)/tankVG, (1.4*tankMGAir+1.329*tankMGH2o)/(tankMGAir+tankMGH2o), 0.85, tankOutletDia); % slightly larger discharge coef to compensate for internal tank contour...
            tankOutletVFlow = tankOutletMFlow/tankLRho;
            mDot = mDot + tankOutletMFlow;
        end
        if(mDot < 0)
            tankMGAir = tankMGAir + mDot*(tankMGAir/(tankMGAir+tankMGH2o))*dt;
            tankMGH2o = tankMGH2o + mDot*(tankMGH2o/(tankMGAir+tankMGH2o))*dt;
        else
            tankMGAir = tankMGAir + mDot*dt;
        end
        %
        if(tankMGAir < 0)
            tankMGAir = 0;
        end
        if(tankMGH2o < 0)
            tankMGH2o = 0;
        end
        if(tankMGAir == 0 && tankMGH2o == 0)
            tankPAir = chmbrPAir;
            tankPH2o = chmbrPNvc;
            tankVG = 0;
            ventVel = 0;
            simData.events.values = [simData.events.values; [simTime, "No gas volume in H2O tank!"]];
        else
            tankPAir = tankMGAir*287.05*tankTG/tankVG;
            tankPH2o = tankMGH2o*461.52*tankTG/tankVG;
            % Thermal effect
            ventVel = mDot / ((tankMGAir+tankMGH2o)/tankVG*tankVentDia);
            dT = 0.5*ventVel^2 * mDot * (1/(h2oCp(tankTL)*tankML));
            dT = dT + (0.04*tankSAL*(tankTL-tankTG)*(1/(h2oCp(tankTL)*tankML)));
            tankTG = tankTG + dT*dt;
        end
    elseif((ctrlTankVent == 0 || ctrlTankRun == 0) && tankVG == 0)
        if(ctrlTankVent == 0)
            mDot = ffH2oVent*mDotThruOrifice(ambientP, tankPAir+tankPH2o, tankLRho, h2oCp(tankTL)/h2oCv(tankTL), 0.8, tankVentDia);
        end
        if(ctrlTankRun == 0)
            tankPBack = backPress(tankPAir+tankPH2o, tankOutletMFlow)-1;
            tankOutletMFlow = ffH2oOutlet*mDotThruOrifice(tankPBack, tankPAir+tankPH2o, (tankML)/tankVL, 1.327, 0.85, tankOutletDia); % slightly larger discharge coef to compensate for internal tank contour...
            tankOutletVFlow = tankOutletMFlow/tankLRho;
            mDot = mDot + tankOutletMFlow;
        end
        if(mDot <= 0)
            tankML = tankML + mDot*dt;
        else
            fprintf("H2O Tank sucking!\n");
        end
        if(tankML < 0)
            tankML = 0;
            simData.events.values = [simData.events.values; [simTime, "All H2O gone!"]];
        else
            tankVL = tankVL + (mDot/tankLRho)*dt;
            tankVG = tankVG - (mDot/tankLRho)*dt;
        end
    end
    %% Diaphragm
    dV = 0;
    %
    if(tankVG > 0 && chmbrVG > 0)
        if(tankVG - abs((mDot/tankLRho)*dt) ~= 0)
            chmbrPAir = chmbrMGAir*287.05*chmbrTG/chmbrVG;
            chmbrPNvc = chmbrMGNvc*33.25*chmbrTG/chmbrVG;
            tankPAir = tankMGAir*287.05*tankTG/tankVG;
            tankPH2o = tankMGH2o*461.52*tankTG/tankVG;
            dV = vDiaph((chmbrPAir+chmbrPNvc)-(tankPAir+tankPH2o), dPLast, KIDiaph, dt, 0.0001, 0.9, 0.02, ffDiaph);
            dPLast = (chmbrPAir+chmbrPNvc)-(tankPAir+tankPH2o);
        else
            tankPAir = chmbrPAir;
            tankPH2o = chmbrPNvc;
            tankVG = 0;
            dV = 0;
        end
    elseif(chmbrVG <= 0)
        chmbrPAir = tankPAir;
        chmbrPNvc = tankPH2o;
        chmbrVG = 0;
        dV = 0;
    elseif(tankVG <= 0)
        tankPAir = chmbrPAir;
        tankPH2o = chmbrPNvc;
        tankVG = 0;
        dV = 0;
    else
        fprintf("Critical Error: No gas volume in chamber or tank\n");
        simData.events.values = [simData.events.values; [simTime, "No gas volume in chamber or tank!"]];
    end
    chmbrVG = chmbrVG + dV*dt;
    tankVG = tankVG - dV*dt;
    if(chmbrVG < 0)
        chmbrVG = 0;
    end
    %
    %% Simulation
    if livePlot
        yyaxis left
        plot(simTime, (chmbrPAir+chmbrPNvc)/1000, '.b', 'MarkerSize', 2);
        plot(simTime, ambientP/1000, '.r', 'MarkerSize', 2);
        plot(simTime, nvcVP(chmbrTL)/1000, '.c', 'MarkerSize', 1);
        yyaxis right
        plot(simTime, chmbrTL, '.g', 'MarkerSize', 2);
        plot(simTime, chmbrTG, '.y', 'MarkerSize', 2);
        drawnow;
    end
    if stoData
        simData.time = [simData.time; simTime];
        
        simData.physics.mass.labels = ["chmbrML", "chmbrMGAir", "chmbrMGNvc", "tankML", "tankMGAir", "tankMGH2o"];
        newMRow = [chmbrML, chmbrMGAir, chmbrMGNvc, tankML, tankMGAir, tankMGH2o];
        
        simData.physics.pressure.labels = ["ambientP/1000", "(chmbrPNvc+chmbrPAir)/1000", "nvcVP(chmbrTL)/1000", "(tankPH2o+tankPAir)/1000", "h2oVP(tankTL)/1000"];
        newPRow = [ambientP/1000, (chmbrPNvc+chmbrPAir)/1000, nvcVP(chmbrTL)/1000, (tankPH2o+tankPAir)/1000, h2oVP(tankTL)/1000];
       
        simData.physics.temperature.labels = ["chmbrTG", "chmbrTL", "tankTG", "tankTL"];
        newTRow = [chmbrTG, chmbrTL, tankTG, tankTL];
        
        simData.physics.volume.labels = ["chmbrVG", "chmbrVL", "tankVG", "tankVL", "diaphDV"];
        newVRow = [chmbrVG, chmbrVL, tankVG, tankVL, dV];
        
        if simTime == 0
            if size(simData.physics.mass.values, 2) ~= size(newMRow, 2)
                simData.physics.mass.values = zeros([0, size(newMRow, 2)]);
            end
            if size(simData.physics.pressure.values, 2) ~= size(newPRow, 2)
                simData.physics.pressure.values = zeros([0, size(newPRow, 2)]);
            end
            if size(simData.physics.temperature.values, 2) ~= size(newTRow, 2)
                simData.physics.temperature.values = zeros([0, size(newTRow, 2)]);
            end
            if size(simData.physics.volume.values, 2) ~= size(newVRow, 2)
                simData.physics.volume.values = zeros([0, size(newVRow, 2)]);
            end
        end
        simData.physics.mass.values = [simData.physics.mass.values; newMRow];
        simData.physics.pressure.values = [simData.physics.pressure.values; newPRow];
        simData.physics.temperature.values = [simData.physics.temperature.values; newTRow];
        simData.physics.volume.values = [simData.physics.volume.values; newVRow];
    end
    fprintf("%0.2f s\n", simTime)
    
end
hold off;

figure(1)
plot(simData.time, simData.physics.pressure.values);
legend(simData.physics.pressure.labels);
xlim([0, simTime]);
%
figure(2)
plot(simData.time, simData.physics.volume.values(:,1:end));
legend(simData.physics.volume.labels(1:end));
xlim([0, simTime]);
figure(3)
plot(simData.time, simData.physics.mass.values(:,4));
yyaxis right;
plot(simData.time, simData.physics.mass.values(:,5)+simData.physics.mass.values(:,6));
yyaxis left;
legend(simData.physics.mass.labels(4), "tankMGTotal");
%ylim([282 300])
xlim([0, simTime]);
%
figure(1)



%H2O Tank Back Pressure function (assumed dependent vars are total tank pressure, flow rate out of tank)
function tankPBack = backPress(tankPTotal, mDot)
    tankPBack = tankPTotal; %assume NO PRESSURE DROP! (only force is diaphragm pushing on water)
end
%Diaphragm dV characteristic (volume defined as internal "bulb" volume, negative value indicates diaphragm has flipped "inside-out"
function dV = vDiaph(dP, dPLast, KIDiaph, dt, KP, KI, KD, ff) % (dp ~ Pa, dV ~ m^3/s)
    KIDiaph = KIDiaph + KI*(dP*dt);
    dV = ff * (KP*dP + KIDiaph - KD*(dP-dPLast)/dt);
end
%Evaporation characteristic of H2O (surfArea in m^3, gasP in Pa, gasT in K, fluid T in K)
function mDot = h2oEvap(sigE, sigC, surfArea, gasP, gasT, fluidT) %a positive (+) mDot signifies mass leaving liquid (evaporating)
    mDot = surfArea * sqrt((2.991507625e-26)/(2*pi*physconst('boltzmann')))*(sigE*h2oVP(fluidT)/sqrt(fluidT) - sigC*gasP/sqrt(gasT));
end
%Vapor Pressure of H2O (T in K, vp in Pa)
function vp = h2oVP(T)
    if(T >= 255.9 && T <= 323)
        A = 4.6543;
        B = 1435.264;
        C = -64.848;
    elseif(T >= 323 && T <= 298)
        A = 5.40221;
        B = 1838.675;
        C = -31.737;
    elseif(T >= 298 && T <= 343.5)
        A = 6.20963;
        B = 2354.731;
        C = 7.559;
    elseif(T >= 343.5 && T <= 373)
        A = 5.08354;
        B = 1663.125;
        C = -45.622;
    elseif(T > 373 && T <= 573)
        A = 3.55959;
        B = 643.748;
        C = -198.043;
    end
    vp = 100000 * 10^(A-(B/(T+C)));
end
%Enthalpy of Vaporization for H2O (T in K, H in kJ/kg)
function H = h2oHvap(T)
    if(T >= 273 && T <= 373.15)
        A = 0;
        B = 0;
        C = 0;
        D = -0.0000128461;
        E = 0.0111289;
        F = -5.58172;
        G = 3457.04;
    elseif(T > 373.15 && T <= 533.15)
        A = 0;
        B = 0;
        C = -0.0000000760531;
        D = 0.000109752;
        E = -0.0636423;
        F = 14.8365;
        G = 1353.78;
    elseif(T > 533.15 && T <= 647)
        A = -0.0000000092542684094;
        B = 0.000031738114185;
        C = -0.045308375811;
        D = 34.461939143;
        E = -14729.409711;
        F = 3354197.6281;
        G = -317931881.15;
    end
    H = A*T^6 + B*T^5 + C*T^4 + D*T^3 + E*T^2 + F*T + G;
end
%Cp Specific Heat for H2O (T in K, Cp in kJ/(kg K))
function Cp = h2oCp(T)
    if(T >= 273 && T <= 373)
        A = -0.00000000000944814;
        B = 0.0000000175084;
        C = -0.0000128895;
        D = 0.00472286;
        E = -0.862214;
        F = 66.9445;
    elseif(T > 373 && T <= 553)
        A = 0;
        B = 0.00000000120945;
        C = -0.00000205829;
        D = 0.00133141;
        E = -0.384920;
        F = 45.9647;
    elseif(T > 553 && T <= 633.15)
        A = 0.0000000085726562494;
        B = -0.000024465214647;
        C = 0.027918109258;
        D = -15.922976704;
        E = 4538.9410095;
        F = -517315.61835;
    end
    Cp = A*T^5 + B*T^4 + C*T^3 + D*T^2 + E*T + F;
    
end

%Cv Specific Heat for H2O (T in K, Cp in kJ/(kg K))
function Cv = h2oCv(T)
    if(T >= 273.16 && T <= 350)
        A = -0.0000000007182627;
        B = 0.000001068461;
        C = -0.0005977988;
        D = 0.1437213;
        E = -8.214461;
    elseif(T > 350 && T <= 500)
        A = -0.0000000001269975;
        B = 0.0000002321628;
        C = -0.0001500751;
        D = 0.03624206;
        E = 1.540572;
    elseif(T > 500 && T <= 633.15)
        A = 0.000000006615767;
        B = -0.00001440192;
        C = 0.01175154;
        D = -4.261921;
        E = 583.0492;
    end
    Cv = A*T^4 + B*T^3 + C*T^2 + D*T + E;        
end

%Evaporation characteristic of NV 7100 (surfArea in m^3, gasP in Pa, gasT in K, fluid T in K)
function mDot = nvcEvap(sigE, sigC, surfArea, gasP, gasT, fluidT) %a positive (+) mDot signifies mass leaving liquid (evaporating)
    mDot = surfArea * sqrt((4.1528e-25)/(2*pi*physconst('boltzmann')))*(sigE*nvcVP(fluidT)/sqrt(fluidT) - sigC*gasP/sqrt(gasT));
end
%Vapor Pressure of NV 7100 (T in K, vp in Pa)
function vp = nvcVP(T)
    vp = exp(22.415 - 3641.9 * (1/T));
end
%Density of NV 7100 Liquid (T in K)
function dD = nvcRho(T)
    dD = 1.5383 - 0.002269*(T-273.15);
    dD = dD / (1000 * 0.000001); %kg/m^3
end
