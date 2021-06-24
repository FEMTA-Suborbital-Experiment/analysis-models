function [Mw, Kc, Cp, k] = AirProperties(T)

	Ru = 8.314;

	xN2  = 78.084 / 100;
	xO2  = 20.947 / 100;
	xAr  = 0.9340 / 100;
	xCO2 = 0.0350 / 100;	

	[MN2, KcN2, CpN2] 	    = GasProperties('N2', T);
	[MO2, KcO2, CpO2] 	    = GasProperties('O2', T);
	[MAr, KcAr, CpAr] 	    = GasProperties('Ar', T);
	[MCO2, KcCO2, CpCO2] 	= GasProperties('CO2', T);

	Cp   = xN2 * CpN2 + xO2 * CpO2 + xAr * CpAr + xCO2 * CpCO2;
	Mw   = xN2 * MN2  + xO2 * MO2  + xAr * MAr  + xCO2 * MCO2;
	
	Ts   = [83.15 123.15 173.15 198.15 223.15 248.15 258.15 263.15 268.15 273.15 278.15 283.15 288.15 293.15 298.15 303.15 313.15 323.15 333.15 353.15 373.15];
	Kcs  = [7.82 11.69 16.2 18.34 20.41 22.41 23.2 23.59 23.97 24.36 24.74 25.12 25.5 25.87 26.24 26.62 27.35 28.08 28.8 30.23 31.62] ./ 1000;

	Kc   = interp1(Ts, Kcs, T);

	Rair = Ru / Mw;
	k    = Cp / (Cp - Rair);
end