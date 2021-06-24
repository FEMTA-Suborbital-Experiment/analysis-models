classdef AlcoholTank < handle

	properties

		dt;
		res;

		%.. geo properties
		D;
		thk;
		hl;
		hg;

		Nx;
		Ny;

		%.. material matrix
		M;

		%.. diffusivity matrix
		A;

		%.. temperature gradient matrix
		dT;

		%.. temperature matrix;
		T;
	end

	methods
	
		function self = AlcoholTank(res, D, thk, hl, hg, Ti, dt)

			self.dt  = dt;
			self.res = res;
			self.D   = D;
			self.thk = thk;
			self.hl  = hl;
			self.hg  = hg;

			Lx = 2 * self.thk + self.D;
			Ly = 2 * self.thk + self.hg + self.hl;

			self.Nx = floor(Lx * self.res) + 3;
			self.Ny = floor(Ly * self.res) + 3;

			self.M = zeros(self.Nx, self.Ny);

			self.M(1 : self.Nx, 1 : self.Ny) = 3;
			self.M(1, :)                     = 4;
			self.M(self.Nx, :)               = 4;
			self.M(:, 1)                     = 4;
			self.M(:, self.Ny)               = 4;

			x = 2;
			y = 2;

			%.. 1 = PEEK, 2 = ALC, 3 = GAS, 4 = AMB
			for lx = [0 : 1 / self.res : Lx]
				
				y = 2;
				
				for ly = [0 : 1 / self.res : Ly]

					if lx >= self.thk & lx <= (self.thk + self.D) & ly >= (self.thk + self.hg) & ly <= (self.thk + self.hg + self.hl)
						self.M(x, y) = 2;
					elseif (lx < self.thk | lx > self.thk + self.D) | (ly < self.thk | ly > self.thk + self.hg + self.hl)
						self.M(x, y) = 1;
					end

					y = y + 1;
				end

				x = x + 1;
			end

			self.T = zeros(self.Nx, self.Ny);

			for dex = [1 : 1 : 4]
				self.T(find(self.M == dex)) = Ti;
			end
			
			self.A = zeros(self.Nx, self.Ny);
			self.dT = zeros(self.Nx, self.Ny);
		end

		function VentUpdate(self, ks, rhos, cps, lstrho, gam)

			for dex = [1 : 1 : 3]
				self.A(find(self.M == dex)) = ks(dex) / (rhos(dex) * cps(dex));
			end

			self.A(find(self.M == 4)) = ks(1) / (rhos(1) * cps(1));

			self.dT = self.A .* self.res .* del2(self.T) .* self.dt;
			self.dT(find(self.M == 3)) = self.dT(find(self.M == 3)) - self.T(find(self.M == 3)) .* (1 - (rhos(3) / lstrho) ^ (gam - 1));

			self.T(2 : self.Nx - 1, 2 : self.Ny - 1) = self.T(2 : self.Nx - 1, 2 : self.Ny - 1) + self.dT(2 : self.Nx - 1, 2 : self.Ny - 1);
		end

		function ZerogUpdate(self, ks, rhos, cps, Hv, mf, dt)

			for dex = [1 : 1 : 3]
				self.A(find(self.M == dex)) = ks(dex) / (rhos(dex) * cps(dex));
			end

			self.A(find(self.M == 4)) = ks(1) / (rhos(1) * cps(1));

			self.dT = self.A .* self.res .* del2(self.T) .* self.dt;
			
			dq = Hv * mf;
			self.dT(find(self.M == 2)) = self.dT(find(self.M == 2)) - dq / (rhos(2) * (pi/4 * self.D^2 * self.hl) * cps(2)) * dt;

			self.T(2 : self.Nx - 1, 2 : self.Ny - 1) = self.T(2 : self.Nx - 1, 2 : self.Ny - 1) + self.dT(2 : self.Nx - 1, 2 : self.Ny - 1);
		end

		function Tavg = AverageTemp(self, sector)

			Tavg = sum(self.T(find(self.M == sector))) / length(self.T(find(self.M == sector)));
		end

		function ShowTempProfile(self, lims)

			imagesc(self.T);
			colorbar;
			caxis(lims);
			colormap jet;
			drawnow;
		end

		function ShowMesh(self)

			imagesc(self.M);
			colormap default;
			drawnow;
		end

		function ShowDiffusivity(self)

			imagesc(self.A);
			colorbar;
			colormap jet;
			drawnow;
		end
	end
end