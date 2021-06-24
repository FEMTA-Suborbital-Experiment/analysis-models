function [MwM, KcM, rhoM, PvM, CpM, HvM] = AlcoholMixture(alcohols, Xs, T)

	MwM  = 0;
	KcM  = 0;
	rhoM = 0;
	PvM  = 0;
	CpM  = 0;
	HvM  = 0;

	for d = [1 : 1 : length(Xs)]

		[Mw, Kc, rho, Pv, Cp, Hv] = LiquidProperties(alcohols(d), T);

		MwM  = MwM + Mw * Xs(d);
		KcM  = KcM + Kc * Xs(d);
		rhoM = rhoM + rho * Xs(d);
		PvM  = PvM + Pv * Xs(d);
		CpM  = CpM + Cp * Xs(d);
		HvM  = HvM + Hv * Xs(d);
	end

	MwM  = MwM / sum(Xs);
	KcM  = KcM / sum(Xs);
	rhoM = rhoM / sum(Xs);
	PvM  = PvM / sum(Xs);
	CpM  = CpM / sum(Xs);
	HvM  = HvM / sum(Xs);
end    
