classdef profiles
	methods(Static)

		function [index_profile, x_grid, y_grid] = get_alpha_law_graded_profile(nx, ny, dx, dy, relative_index_diff_delta, a_core_radius_m, n0_core, alpha_power)
			% Source: Gloge, D., & Marcatili, E. A. J. (1973). Multimode Theory of Graded‐Core Fibers. Bell System Technical Journal, 52(9), 1563–1578.

			x_grid = linspace(-nx*dx/2, nx*dx/2, nx);
			y_grid = linspace(-ny*dy/2, ny*dy/2, ny);

			n_cladding = n0_core * sqrt(1 - 2*relative_index_diff_delta);
			index_profile = n_cladding * ones(nx, ny);

			for xx = 1:nx
				for yy = 1:ny
					ro = sqrt(x_grid(xx)^2 + y_grid(yy)^2);
					if ro <= a_core_radius_m
						% inside the core
						index_profile(yy,xx) = n0_core * sqrt(1 - 2*relative_index_diff_delta*(ro/a_core_radius_m)^alpha_power);
					else
						% inside the cladding
						index_profile(yy,xx) = n_cladding;
					end
				end
			end
		end

		function [index_rho, rho_arr] = get_alpha_law_graded_profile_radial(nr, dr, relative_index_diff_delta, a_core_radius_m, n0_core, alpha_power)
			% Source: Gloge, D., & Marcatili, E. A. J. (1973). Multimode Theory of Graded‐Core Fibers. Bell System Technical Journal, 52(9), 1563–1578.
			% nr = number of radial points to sample
			% dr = distance between consecutive radial points

			rho_arr = linspace(0,nr*dr,nr);

			n_cladding = n0_core * sqrt(1 - 2*relative_index_diff_delta);
			index_rho = n_cladding * ones(1, nr);

			for ii = 1:nr
				rho = rho_arr(ii);
				if rho <= a_core_radius_m
					% inside the core
					index_rho(ii) = n0_core * sqrt(1 - 2*relative_index_diff_delta*(rho/a_core_radius_m)^alpha_power);
				else
					% inside the cladding
					index_rho(ii) = n_cladding;
				end
			end
		end

		function [index_profile, x_grid, y_grid] = get_step_index_profile(nx, ny, dx, dy, relative_index_diff_delta, a_core_radius_m, n0_core)
			% Source: Kahn EE247 Lecture 2: Geometrical Optics Description of Optical Fibers.
			x_grid = linspace(-nx*dx/2, nx*dx/2, nx);
			y_grid = linspace(-ny*dy/2, ny*dy/2, ny);

			n_cladding = n0_core * (1 - relative_index_diff_delta); %TODO: why is there a factor of 2 difference?
			index_profile = n_cladding * ones(nx, ny);

			for xx = 1:nx
				for yy = 1:ny
					ro = sqrt(x_grid(xx)^2 + y_grid(yy)^2);
					if ro <= a_core_radius_m
						% inside the core
						index_profile(yy,xx) = n0_core;
					else
						% inside the cladding
						index_profile(yy,xx) = n_cladding;
					end
				end
			end
		end

		function [index_rho, rho_arr] = get_step_index_profile_radial(nr, dr, relative_index_diff_delta, a_core_radius_m, n0_core)
			% Source: Kahn EE247 Lecture 2: Geometrical Optics Description of Optical Fibers.
			% nr = number of radial points to sample
			% dr = distance between consecutive radial points

			rho_arr = linspace(0,nr*dr,nr);

			n_cladding = n0_core * (1 - relative_index_diff_delta);
			index_rho = n_cladding * ones(1,nr);

		
			for ii = 1:nr
				rho = rho_arr(ii);
				if rho <= a_core_radius_m
					% inside the core
					index_rho(ii) = n0_core;
				else
					% inside the cladding
					index_rho(ii) = n_cladding;
				end
			end
		end

		function [index_profile, x_grid, y_grid] = get_GIGDC_profile(nx, ny, dx, dy, relative_index_diff_delta, a_core_radius_m, n0_core, alpha_power, dip)
			% Source: Ishigure, T., Endo, H., Ohdoko, K., & Koike, Y. (2004). High-bandwidth plastic optical fiber with W-refractive index profile. IEEE Photonics Technology Letters, 16(9), 2081–2083. 

			x_grid = linspace(-nx*dx/2, nx*dx/2, nx);
			y_grid = linspace(-ny*dy/2, ny*dy/2, ny);

			n_cladding = n0_core * sqrt(1 - 2*relative_index_diff_delta);
			index_profile = n_cladding * ones(nx, ny);

			for xx = 1:nx
				for yy = 1:ny
					ro = sqrt(x_grid(xx)^2 + y_grid(yy)^2);
					if ro <= a_core_radius_m
						% inside the core
						index_profile(yy,xx) = n0_core * sqrt(1 - 2*relative_index_diff_delta*dip*(ro/a_core_radius_m)^alpha_power);
					else
						% inside the cladding
						index_profile(yy,xx) = n_cladding;
					end
				end
			end
		end

		function [index_profile, x_grid, y_grid] = get_multicore_profile(nx, ny, dx, dy, n_cores, core_radii_m, core_xlocs, core_ylocs, n0_cores, relative_index_diff_delta)
			x_grid = linspace(-nx*dx/2, nx*dx/2, nx);
			y_grid = linspace(-ny*dy/2, ny*dy/2, ny);
			n_cladding = mean(n0_cores) * (1 - relative_index_diff_delta); % TODO: is there a factor of 2 for deltan for step index fibers?? is the mean() correct if there are heterogeneous cores?
			index_profile = n_cladding * ones(nx, ny);

			for nn = 1:n_cores
				for xx = 1:nx
					for yy = 1:ny
						ro = sqrt((x_grid(xx)-core_xlocs(nn))^2 + (y_grid(yy)-core_ylocs(nn))^2);
						if ro <= core_radii_m(nn)
							% inside the core
							index_profile(yy,xx) = n0_cores(nn);
						end
					end
				end
			end
		end

	end


end