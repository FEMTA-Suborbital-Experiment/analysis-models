function [Mw, Kc, rho, Pv, Cp, Hv] = LiquidProperties(liquid, T)

    Ru    = 8.314;
    T_ref = 300;
    T_avg = [200 : 1 : 400];

	%.. liquid methanol
	if strcmp(liquid, 'methanol')   

        Mw  = 32.04 / 1000;                           % molar mass           [kg/mol]
        Kc  = 0.203;                                  % thermal conductivity [W/m-k]
        rho = 792.0;                                  % density              [kg/m^3]
        Cp  = 80.0 / Mw;                              % specific heat        [J/kg-k]

        hv_coeff = [45.3, -0.31, 0.4241, 512.6];

        if T <= 365
            pv_coeff = [5.20409 1581.341 -33.50];
        elseif T > 365 && T <= 513.91
            pv_coeff =[5.31301 1676.569 -21.728];
        else 
            pv_coeff = [NaN, NaN, NaN];
        end

    elseif strcmp(liquid, 'ethanol') 

    	Mw  = 46.07 / 1000;                           % molar mass           [kg/mol]
        Kc  = 0.171;                                  % thermal conductivity [W/m-k]
        rho = 789.0;                                  % density              [kg/m^3]
        Cp  = 110.0 / Mw;                             % specific heat        [J/kg-k]

        hv_coeff = [50.43, -0.4475, 0.4989, 513.9];

        if T <= 350
            pv_coeff = [5.37229, 1670.409, -40.191];
        elseif T > 350 && T <= 513.91
            pv_coeff =[4.92531, 1432.526, -61.819];
        else 
            pv_coeff = [NaN, NaN, NaN];
        end

    elseif strcmp(liquid, 'propanol')

    	Mw  = 60.09 / 1000;                           % molar mass           [kg/mol]
        Kc  = 0.154;                                  % thermal conductivity [W/m-k]
        rho = 803.0;                                  % density              [kg/m^3]
        Cp  = 140.0 / Mw;                             % specific heat        [J/kg-k]

        hv_coeff = [52.06, -0.8386, 0.6888, 536.7];

        if T <= 370
            pv_coeff = [5.31384, 1690.864, -51.804];
        elseif T > 370 && T <= 536.71
            pv_coeff = [4.59871, 1300.491, -86.3641];
        else 
            pv_coeff = [NaN, NaN, NaN];
        end

    elseif strcmp(liquid, 'isopropyl')

    	Mw  = 60.10 / 1000;                           % molar mass           [kg/mol]
        Kc  = 0.140;                                  % thermal conductivity [W/m-k]
        rho = 786.0;                                  % density              [kg/m^3]
        Cp  = 3680.0;                                 % specific heat        [J/kg-k]

        hv_coeff = [53.38, -0.708, 0.6538, 508.3];

        if T <= 500
            pv_coeff = [4.57795, 1221.423, -87.474];
        else 
            pv_coeff = [NaN, NaN, NaN];
        end

    elseif strcmp(liquid, 'butanol')

    	Mw  = 74.10 / 1000;                           % molar mass           [kg/mol]
        Kc  = 0.154;                                  % thermal conductivity [W/m-k]
        rho = 810.0;                                  % density              [kg/m^3]
        Cp  = 175.0 / Mw;                             % specific heat        [J/kg-k]

        hv_coeff = [62.53, -0.6584, 0.696, 562.9];

        if T <= 390
            pv_coeff = [4.54607, 1351.555, -93.34];
        elseif T > 390 && T <= 562.98
            pv_coeff = [4.42921, 1305.001, -94.676];
        else 
            pv_coeff = [NaN, NaN, NaN];
        end
        
    elseif strcmp(liquid, 'pentanol')

    	Mw  = 88.15 / 1000;                           % molar mass           [kg/mol]
        Kc  = 0.153;                                  % thermal conductivity [W/m-k]
        rho = 814.0;                                  % density              [kg/m^3]
        Cp  = 205.0 / Mw;                             % specific heat        [J/kg-k]

        hv_coeff = [67.55, -0.8195, 0.8272, 588.2];

        if T <= 410
            pv_coeff = [4.68277, 1492.549, -91.621];
        elseif T > 410 & T <= 513.79
            pv_coeff = [3.97383, 1106.11, -134.578];
        else 
            pv_coeff = [NaN, NaN, NaN];
        end
    end

    %.. heat of vaporization [J/kg]
    Tr = T / hv_coeff(4);                                  
    Hv = hv_coeff(1) * exp(-hv_coeff(2) * Tr) * (1 - Tr) ^ hv_coeff(3);               
    Hv = Hv * 1000 / Mw;                                        

    %.. vapor pressure [pa]
    Pv = 10 ^ (pv_coeff(1) - pv_coeff(2) / (T + pv_coeff(3))) * 10 ^ 5;
end