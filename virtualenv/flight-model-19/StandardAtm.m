    
function [T, p, rho] = StandardAtm(h)
    
    rho_sea = 1.2250;     % density at sea level, kg/ m^3
    p_sea   = 1.01325E5;  % pressure at sea level, N/m^2
    R       = 287;        % gas constant, J/kg/K
    g_zero  = 9.81;       % gravitatoinal acceleration, m/s^2

    T_set   = [288.16,216.66,216.66,282.66,282.66,165.66,165.66,225.66]; % list of temperature points that define endpoints of each layer (starting at the ground), K
    h_set   = [0,11000,25000,47000,53000,79000,90000,105000]; %list of altitude points that define endpoints of each layer (starting at the ground), m
    a_set   = [-6.5*10^-3,0,3*10^-3,0,-4.5*10^-3,0,4*10^-3];%list of gradient layer slopes (starting at the ground), K/m
    p_set   = [p_sea, 2.27e4, 2.5273e3, 1.2558e2, 61.493, 1.2595, 0.162723]; % list of pressure at each layer endpoint, N/m^2
    rho_set = [rho_sea, 3.648e-1, 4.0639e-2, 1.5535e-3, 7.5791e-4, 2.184e-5, 1.84114e-6]; % list of density at each layer endpoint, kg/m^3
            
    if h <= h_set(2)
        layer = 1;
    elseif h <= h_set(3)
        layer = 2;
    elseif h <= h_set(4)
        layer = 3;
    elseif h <= h_set(5)
        layer = 4;
    elseif h <= h_set(6)
        layer = 5;
    elseif h <= h_set(7)
        layer = 6;
    %elseif h <= h_set(8)
    else
        layer = 7;
    end
    
    
    if mod(layer, 2) == 1
        %Gradient layer
        T = T_set(layer) + (a_set(layer) * (h - h_set(layer))); % temperature equation for gradient layer, K
        p = p_set(layer) * ((T / T_set(layer)) ^ (-g_zero / (a_set(layer) * R))); % pressure equation for gradient layer, N/m^2
        rho = rho_set(layer)* ((T / T_set(layer))^((-g_zero / (a_set(layer)*R))+1)); % density equation for gradient layer, kg/m^3
        
    else
        %Isothermal layers
        T = T_set(layer); % temperature for isothermal layer, K
        p = p_set(layer) * exp((-g_zero * (h - h_set(layer))) / (R * T)); % pressure equation for isothermal layer, N/m^2
        rho = (p * rho_set(layer)) / p_set(layer); % density equation for isothermal layer, kg/m^3
    end
    
end
       