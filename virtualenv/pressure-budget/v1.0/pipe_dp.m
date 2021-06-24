function [out1] = pipe_dp(in1,in2,in3,in4,in5,in6,in7,in8)
    %Pipe delta-P
    %   Finds deltaP in any pipe length with optional given bend radius
    mF = in7; %moody friction factor
    rho = in5; %kg/m^3 fluid density
    mu = in6; %m/s (fluid flow velocity)
    bR = in2; %m bend radius
    theta = in1; %bend angle (deg)
    fCoef = in8; %pipe friction coefficient
    L = in3; %m pipe length
    iD = in4; %m pipe inside diameter
    if(theta ~= 0 && bR ~= 0)
        % bendLossCoef.csv is the look up table for the bend loss coefficient. More at the URL provided below
        lossCoefTable = readtable('bendLossCoef.csv'); %bend loss coefficient -> http://www.thermopedia.com/content/5000/80B(FAPDI)Fig3.gif
        lossCoefTable = lossCoefTable{:,:};
        xVal = bR/iD;
        if xVal > 10
            xVal = 10;
        elseif xVal < 0.5
            fprintf("Out of Bounds Loss Coefficient\n");
        end
        %following block allows for exterpolation
        upperLim = ceil(xVal * 10)/10;
        lowerLim = floor(xVal * 10)/10;
        upperCol = find(lossCoefTable(1,:)*10==int8(upperLim*10));
        while isempty(upperCol) && upperLim < 10
            upperLim = upperLim + 0.1;
            upperCol = find(lossCoefTable(1,:)*10==int8(upperLim*10));
        end
        lowerCol = find(lossCoefTable(1,:)*10==int8(lowerLim*10));
        while isempty(lowerCol) && lowerLim < 10
            lowerLim = lowerLim - 0.1;
            lowerCol = find(lossCoefTable(1,:)*10==int8(lowerLim*10));
        end
        %
        upperLim = ceil(theta);
        upperLim = ceil(theta/10)*10;
        lowerLim = floor(theta);
        lowerLim = floor(theta/10)*10;
        upperRow = find(lossCoefTable(:,1)==int16(upperLim));
        while isempty(upperRow) && upperLim < 190
            upperLim = upperLim + 10;
            upperRow = find(lossCoefTable(:,1)==int16(upperLim));
        end
        lowerRow = find(lossCoefTable(:,1)==int16(lowerLim));
        while isempty(lowerRow) && lowerLim < 190
            lowerLim = lowerLim + 10;
            lowerRow = find(lossCoefTable(:,1)==int16(lowerLim));
        end
        %
        kb = (lossCoefTable(upperRow,upperCol) + lossCoefTable(lowerRow,lowerCol)) / 2;
        %dP for bent pipe
        dP = (0.5 * mF * rho * (mu*mu) * ((pi*bR)/iD) * (theta/180)) + (0.5 * kb * rho * (mu*mu));
    else
        %dp for straight pipe
        dP = fCoef * (L/iD) * (rho/2) * (mu*mu);
    end
    out1 = dP;
end

