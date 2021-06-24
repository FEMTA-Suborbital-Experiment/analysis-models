function [mDot] = mDotThruOrifice(in1,in2,in3,in4,in5,in6)
%mDotThruOrifice calculates mass flow (mDot) in kg/s through an orifice
%given pressures and working fluid properties
%
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
