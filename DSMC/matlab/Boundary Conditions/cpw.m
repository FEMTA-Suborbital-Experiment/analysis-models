function [Cp] = cpw(T)

	Ts = [175 200 225 250 275 300 325 350 375 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200 1250 1300 1350 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3500 4000 4500 5000 5500 6000];
	Cps = [1.85 1.851 1.852 1.855 1.859 1.864 1.871 1.88 1.89 1.901 1.926 1.954 1.984 2.015 2.047 2.08 2.113 2.147 2.182 2.217 2.252 2.288 2.323 2.358 2.392 2.425 2.458 2.49 2.521 2.552 2.609 2.662 2.711 2.756 2.798 2.836 2.872 2.904 2.934 2.962 2.987 3.011 3.033 3.053 3.072 3.09 3.163 3.217 3.258 3.292 3.322 3.35];
    
	if (T < min(Ts))
		Cp = Cps(1) .* 1000;
		return;
    end
    
	Cp = interp1(Ts, Cps, T) .* 1000;
end