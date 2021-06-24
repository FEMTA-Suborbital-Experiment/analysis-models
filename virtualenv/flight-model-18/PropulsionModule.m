classdef PropulsionModule < handle

	properties

		% FIXED PROPERTIES

		Cdv;                             % vent orifice discharge coefficient           [...]
		Cdp;                             % propellant orifice discharge coefficient     [...]

		Av;                              % vent orifice area                            [m^2]
		Ap;                              % propellant orifice area                      [m^2]

		Dac;                             % diameter of the alcohol chamber              [m]
		tac;                             % thermal barrier thickness of alcohol chamber [m]
		hac;                             % total alcohol chamber height                 [m]
		hal;                             % height of alcohol liquid in chamber          [m]

		mal0;                            % initial mass of alcohol liquid               [kg]
		Val0;                            % initial volume of alcohol liquid             [m^3]
		Vag0;                            % initial volume of the alcohol gas            [m^3]
		mp0;                             % initial mass of propellant                   [kg]

		alcohols;                        % array of alcohol names for the mixture       [...]
		xs;                              % array of alcohol mole fractions              [...]

		rhoPK = 1320;                    % density of PEEK plastic                      [kg/m^3]
		KcPK  = 0.201;                   % thermal conductivity of PEEK plastic         [W/m-k]
		CpPK  = 320.0;                   % specific heat capacity of PEEK plastic       [J/kg-k]

		rhop0 = 997.0;                   % initial density of propellant                [kg/m^3]

		MeshRes;                         % mesh resolution                              [nodes/mm]

		% FLIGHT PROFILE

		t = [];                          % time profile                                 [s]
		h = [];                          % altitude profile                             [m]
		a = [];                          % acceleration profile                         [a]

		pThresh;                         % pressure threshold

		dt;                              % time incremenet                              [s]

		% DYNAMIC PROPERTIES

		%.. ambient conditions
		T_am   = [];
		P_am   = [];
		rho_am = [];

		%.. gas chamber properties
		T_ag   = [];
		P_ag   = [];
		rho_ag = [];
		m_ag   = [];
		mf_ag  = [];
		Cp_ag  = [];
		k_ag   = [];

		%.. alcohol liquid properties
		Cp_al  = [];
		Hv_al  = [];
		T_al   = [];

		%.. propellant properties
		P_p    = [];
		m_p    = [];
		mf_p   = [];

		% OBJECTS
		alctank;                         % instance of alcohol tank class               [AlcoholTank]

		% UNIVERSAL CONTANTS

		Rair = 287.058;                  % specific gas constant of air                 [J/kg-k]

	end

	methods

		% flightProfile: [t, h, a]       (s, m, m/s^2)
		% ventPort:      [Cdv, Av]       (.., m^2)
		% propPort:      [Cdp, Ap]       (.., m^2)
		% dims:          [Dac, tac, hac] (m)
		% fluidQtys:     [mp0, mal0]     (g)

		function self = PropulsionModule(flightProfile, pThresh, ...
			                             ventPort, propPort, ...
			                             dims, fluidQtys, ...
			                             alcohols, xs, ...
			                             meshResolution)

			%.. loading flight profile
			self.t  = flightProfile(1, :);
			self.h  = flightProfile(2, :);
			self.a  = flightProfile(3, :);

			self.dt = self.t(2) - self.t(1);

			self.pThresh = pThresh;

			%.. setting port properties
			self.Cdv = ventPort(1);
			self.Av  = ventPort(2);

			self.Cdp = propPort(1);
			self.Ap  = propPort(2);

			%.. setting module dim properties
			self.Dac  = dims(1);
			self.tac  = dims(2);
			self.hac  = dims(3);

			%.. setting module fluid quantities
			self.mp0  = fluidQtys(1) / 1000;
			self.mal0 = fluidQtys(2) / 1000;

			%.. setting alcohol mixture properties
			self.alcohols = alcohols;
			self.xs       = xs;

			self.MeshRes  = meshResolution;
		end

		function Fly(self, outputVideo, thermalRange)

			ReInit(self);

			%.. initial properties
			[self.T_am(1), self.P_am(1), self.rho_am(1)] = StandardAtm(self.h(1));

			self.rho_ag(1) = self.rho_am(1);
			self.P_p(1)    = self.P_am(1);
			self.P_ag(1)   = self.P_am(1);
			self.T_ag(1)   = self.T_am(1);
			self.T_al(1)   = self.T_am(1);
			self.m_ag(1)   = self.Vag0 * self.rho_am(1);
			self.mf_ag(1)  = 0;

			[~, ~, self.Cp_ag(1), self.k_ag(1)] = AirProperties(self.T_am(1));

			[~, ~, ~, ~, self.Cp_al(1), self.Hv_al(1)] = AlcoholMixture(self.alcohols, self.xs, self.T_am(1));

			self.m_p(1)  = self.mp0;
			self.mf_p(1) = 0.0;

			%.. looping

			if outputVideo
				
				fig = figure('visible','off');

				writer = VideoWriter('thermal.mp4','MPEG-4');
				writer.Quality = 100;

				open(writer);
			end

			for n = 2 : length(self.t)

				[self.T_am(n), self.P_am(n), self.rho_am(n)] = StandardAtm(self.h(n));

				if self.P_am(n) >= self.pThresh & self.a(n) ~= 0

					% VENT PHASE

					%.. constant values
					self.m_p(n)  = self.m_p(n - 1);
					self.mf_p(n) = 0;

					%.. air mach number
					if self.P_am(n) < self.P_ag(n - 1)

						% ambient pressure is lower
						P0 = self.P_ag(n - 1);
						P1 = self.P_am(n);
						T0 = self.T_ag(n - 1);

						% mach number of flow
						M = sqrt(2 / (self.k_ag(n - 1)) * ((P1 / P0) ^ (-(self.k_ag(n - 1) - 1) / self.k_ag(n - 1)) - 1));

						if M > 1
							M = 1;
						end

						self.mf_ag(n) = self.Cdv * self.Av * P0 / sqrt(T0) * ...
						                sqrt(self.k_ag(n - 1) / self.Rair) * (1 + (self.k_ag(n - 1) - 1) / 2 * M ^ 2) ^ ...
						                (-(self.k_ag(n - 1) + 1) / (2 * (self.k_ag(n - 1) - 1))) * M;
					else 

						% no induced mass flow rate
						self.mf_ag(n) = 0;
					end

					%.. air mass and density
					dm = self.mf_ag(n) * self.dt;

					rho_min = self.rho_am(1) * (self.P_am(n) / self.P_ag(1)) ^ (1 / self.k_ag(n - 1));
					dm_min  = self.m_ag(n - 1) - rho_min * self.Vag0;

					if dm > dm_min
						dm = dm_min;
					end

					self.m_ag(n)   = self.m_ag(n - 1) - dm;
					self.rho_ag(n) = self.rho_ag(n - 1) * self.m_ag(n) / self.m_ag(n - 1);

					%.. air and propellant pressure
					self.P_ag(n) = self.P_ag(1) * (self.rho_ag(n) / self.rho_ag(1)) ^ (self.k_ag(n - 1));
					self.P_p(n)  = self.P_ag(n);

                    %.. air thermodynamic properties
                    [~, Kc_ag, self.Cp_ag(n), self.k_ag(n)] = AirProperties(self.T_ag(n - 1));

                    %.. alcohol liquid thermodynamic properties
                    [~, Kc_al, rho_al, ~, self.Cp_al(n), self.Hv_al(n)] = AlcoholMixture(self.alcohols, self.xs, self.T_al(n - 1));

                    %.. solve FEM for vent phase
                    VentUpdate(self.alctank, ...
                    	       [self.KcPK,         Kc_al,         Kc_ag], ...
                    	       [self.rhoPK,        rho_al,        self.rho_ag(n)], ...
                    	       [self.CpPK,         self.Cp_al(n), self.Cp_ag(n)], ...
                    	       self.rho_ag(n - 1), self.k_ag(n));

                    %.. updating alc liquid and air temp
                    self.T_al(n) = AverageTemp(self.alctank, 2);
                    self.T_ag(n) = AverageTemp(self.alctank, 3);

				elseif self.P_am(n) < self.pThresh & self.a(n) ~= 0

					% COAST PHASE

					%.. gas chamber properties
					self.rho_ag(n) = self.rho_ag(n - 1);
					self.m_ag(n)   = self.m_ag(n - 1);
					self.mf_ag(n)  = 0;

					%.. propellant properties
					self.m_p(n)  = self.m_p(n - 1);
					self.mf_p(n) = 0;

                    %.. air thermodynamic properties
                    [~, Kc_ag, self.Cp_ag(n), self.k_ag(n)] = AirProperties(self.T_ag(n - 1));

                    %.. alcohol liquid thermodynamic properties
                    [~, Kc_al, rho_al, self.P_ag(n), self.Cp_al(n), self.Hv_al(n)] = AlcoholMixture(self.alcohols, self.xs, self.T_al(n - 1));

                    self.P_p(n) = self.P_ag(n);

                    %.. solve FEM for vent phase
                    VentUpdate(self.alctank, ...
                    	       [self.KcPK,     Kc_al,         Kc_ag], ...
                    	       [self.rhoPK,    rho_al,        self.rho_ag(n)], ...
                    	       [self.CpPK,     self.Cp_al(n), self.Cp_ag(n)], ...
                    	       self.rho_ag(n), self.k_ag(n));

                    %.. updating alc liquid and air temp
                    self.T_al(n) = AverageTemp(self.alctank, 2);
                    self.T_ag(n) = AverageTemp(self.alctank, 3);

				elseif self.P_am(n) < self.pThresh & self.a(n) == 0

					% ZERO-G PHASE

                    %.. air thermodynamic properties
                    [~, Kc_ag, self.Cp_ag(n), self.k_ag(n)] = AirProperties(self.T_ag(n - 1));

                    %.. alcohol liquid thermodynamic properties
                    [~, Kc_al, rho_al, self.P_ag(n), self.Cp_al(n), self.Hv_al(n)] = AlcoholMixture(self.alcohols, self.xs, self.T_al(n - 1));

                    %.. mass flow rate and mass of propellant
                    self.P_p(n)    = self.P_ag(n);
                   	self.mf_p(n)   = 2 * self.Cdp * (2 * self.Ap) * sqrt(2 * self.rhop0 * (self.P_p(n) - self.P_am(n)));

                   	if self.mf_p(n) * self.dt > self.m_p(n - 1)

                   		self.mf_p(n) = self.m_p(n - 1) / self.dt;
					end

                   	self.m_p(n) = self.m_p(n - 1) - self.mf_p(n) * self.dt;
                   	
                   	%.. alcohol gasses
                   	Vf             = self.mf_p(n) / self.rhop0;
                   	self.mf_ag(n)  = Vf * self.rho_ag(n - 1);
                   	self.m_ag(n)   = self.m_ag(n - 1) + self.mf_ag(n) * self.dt;
                   	self.Vag0      = self.Vag0 + Vf * self.dt;
                   	self.rho_ag(n) = self.m_ag(n) / self.Vag0;

                    %.. solve FEM for vent phase
                    ZerogUpdate(self.alctank, ...
                    	        [self.KcPK,    Kc_al,         Kc_ag], ...
                    	        [self.rhoPK,   rho_al,        self.rho_ag(n)], ...
                    	        [self.CpPK,    self.Cp_al(n), self.Cp_ag(n)], ...
                    	        self.Hv_al(n), self.mf_ag(n), self.dt);

                    %.. updating alc liquid and air temp
                    self.T_al(n) = AverageTemp(self.alctank, 2);
                    self.T_ag(n) = AverageTemp(self.alctank, 3);
				end

				%.. writing to video
				if outputVideo

					if mod(n, 100) == 0
						
						fig; imagesc(self.alctank.T);
						fig; colorbar;
						fig; caxis([self.T_am(1) - thermalRange, self.T_am(1)]);
						fig; colormap jet;
						fig; drawnow;

						writeVideo(writer, getframe(fig));

						clc;
						fprintf('frame saved: %.2f s\n', self.t(n));
					end
				end
			end

			if outputVideo

				close(writer);
				close all;
			end
		end

		function ReInit(self)

			%.. ambient conditions
			self.T_am   = [];
			self.P_am   = [];
			self.rho_am = [];

			%.. gas chamber properties
			self.T_ag   = [];
			self.P_ag   = [];
			self.rho_ag = [];
			self.m_ag   = [];
			self.mf_ag  = [];
			self.Cp_ag  = [];
			self.k_ag   = [];

			%.. alcohol liquid properties
			self.Cp_al  = [];
			self.Hv_al  = [];
			self.T_al   = [];

			%.. propellant properties
			self.P_p    = [];
			self.m_p    = [];
			self.mf_p   = [];

			%.. initializing calculated properties
			[T0, ~, ~] = StandardAtm(self.h(1));
			[~, ~, rho_al, ~, ~, ~] = AlcoholMixture(self.alcohols, self.xs, T0);

			self.Val0  = self.mal0 / rho_al;
			self.hal  = (4 * self.Val0) / (pi * self.Dac ^ 2);
			
			if self.hal >= self.hac

				error('error: too much alcohol liquid for tank size');
			end

			self.Vag0 = pi / 4 * self.Dac ^ 2 * (self.hac - self.hal);

			%.. initializing tank FEM
			self.alctank = AlcoholTank(self.MeshRes, ...
				                       self.Dac, ...
				                       self.tac, ...
				                       self.hal, ...
				                       self.hac - self.hal, ...
				                       T0, ...
				                       self.dt);
		end
	end
end