classdef profiles
	methods(Static)

		function [dn_fiber, n_clad, rho_arr, alpha_part, non_alpha_part] = get_nr_MJL_profile(B, nr, dr, relative_dn, a_core_radius_m, n_clad, alpha0)
			% Ming Jun Li's patent. B is coefficients of perturbation in polynomial basis.
			% dn_fiber = index variation in (x,y) that is an offset from the cladding index, i.e. n(x,y) = dn_fiber(x,y) + n_clad
			rho_arr = linspace(0, nr*dr, nr);

			N = length(B); 	% number of correction terms in polynomial basis

            % form generalized perturbation H using correction terms
            G = zeros(1, nr);

            pure_alpha 	= zeros(1, nr);	% pure alpha contributions
            H 			= zeros(1, nr); % non alpha contributions
			for rr = 1:nr
				rho_norm = rho_arr(rr) / a_core_radius_m;
				if rho_norm < 1 % region where correction is valid
					
					% form pure alpha law core
					pure_alpha(rr) = rho_norm^alpha0;

					% add non-alpha contributions
					for nn = 1:N 
						H(rr) = H(rr) + B(nn) * utils.poly_basis(nn, rho_norm);
					end
				end
			end
			G = pure_alpha + H;
			alpha_part = pure_alpha;
			non_alpha_part = H;


			n_core = n_clad / sqrt(1 - 2*relative_dn);
			n_rho = n_clad * ones(1, nr);
			for rr = 1:nr
				ro = rho_arr(rr);
				if ro <= a_core_radius_m
					% inside the core
					n_rho(rr) = n_core * sqrt(1 - 2*relative_dn*G(rr));
				else
					% inside the cladding
					n_rho(rr) = n_clad;
				end
			end
			dn_fiber = n_rho - n_clad;

		end

		function [index_profile, x_arr, y_arr] = synthesize_nxy_from_nr(n_rho, rho_arr, nan_fill_val)
			% rho is the radial variable
			% n_rho is the radial distribution of the index evaluated at rho_arr
			% n_rho and rho_arr are arrays of same length
			% nan_fill_val is the index to put into NaN locations of synthesized index profile

			n_theta = 100;
			theta_arr = linspace(0, 2*pi, n_theta);
		    [R, THETA] = meshgrid(rho_arr, theta_arr);
		    N_RHO = repmat(n_rho, length(theta_arr), 1);

		    x_arr = linspace(-max(rho_arr), max(rho_arr), length(rho_arr));
		    y_arr = x_arr;
		
		    [Xq, Yq] = meshgrid(x_arr, y_arr);
		   
		    index_profile = griddata(R.*cos(THETA), R.*sin(THETA), N_RHO, Xq, Yq);
		    % turn off annoying warning messages
		    w = warning('query','last');
		    id = w.identifier;
		    warning('off',id);

		    % griddata gives nan for extrapolated points so replace NaNs with a refractive index value
		    index_profile(isnan(index_profile)) = nan_fill_val; 
		end



		function [dn_fiber, n_clad, x_grid, y_grid] = get_nxy_alpha_law_graded_profile(nx, ny, dx, dy, relative_dn, a_core_radius_m, n_clad, alpha_power)
			% Source: Physical modeling of 10 gbe optical communication systems - Gholami, Molin and Sillard, JLT 2011
			% dn_fiber = index variation in (x,y) that is an offset from the cladding index, i.e. n(x,y) = dn_fiber(x,y) + n_clad
			x_grid = linspace(-nx*dx/2, nx*dx/2, nx);
			y_grid = linspace(-ny*dy/2, ny*dy/2, ny);

			n_core = n_clad / sqrt(1 - 2*relative_dn);
			n_xy = n_clad * ones(nx, ny);
			for xx = 1:nx
				for yy = 1:ny
					ro = sqrt(x_grid(xx)^2 + y_grid(yy)^2);
					if ro <= a_core_radius_m
						% inside the core
						n_xy(yy,xx) = n_core * sqrt(1 - 2*relative_dn*(ro/a_core_radius_m)^alpha_power);
					else
						% inside the cladding
						n_xy(yy,xx) = n_clad;
					end
				end
			end
			dn_fiber = n_xy - n_clad;
		end

		function [dn_fiber, n_clad, x_grid, y_grid] = get_nxy_multicore_profile(nx, ny, dx, dy, n_cores, core_radii_m, core_xlocs, core_ylocs, n_clad, relative_dn)
			x_grid = linspace(-nx*dx/2, nx*dx/2, nx);
			y_grid = linspace(-ny*dy/2, ny*dy/2, ny);
			
			n_core = n_clad / sqrt(1 - 2*relative_dn);
			n_xy = n_clad * ones(nx, ny);

			for nn = 1:n_cores
				for xx = 1:nx
					for yy = 1:ny
						ro = sqrt((x_grid(xx)-core_xlocs(nn))^2 + (y_grid(yy)-core_ylocs(nn))^2);
						if ro <= core_radii_m(nn)
							% inside the core
							n_xy(yy,xx) = n_core;
						end
					end
				end
			end
			dn_fiber = n_xy - n_clad;
		end

		function [dn_fiber, n_clad, rho_arr] = get_nr_alpha_law_graded_profile(nr, dr, relative_dn, a_core_radius_m, n_clad, alpha_power)
			% Source: Physical modeling of 10 gbe optical communication systems - Gholami, Molin and Sillard, JLT 2011
			% dn_fiber = index variation in (x,y) that is an offset from the cladding index, i.e. n(x,y) = dn_fiber(x,y) + n_clad
			rho_arr = linspace(0, nr*dr, nr);

			n_core = n_clad / sqrt(1 - 2*relative_dn);
			n_rho = n_clad * ones(1, nr);
			for rr = 1:nr
				ro = rho_arr(rr);
				if ro <= a_core_radius_m
					% inside the core
					n_rho(rr) = n_core * sqrt(1 - 2*relative_dn*(ro/a_core_radius_m)^alpha_power);
				else
					% inside the cladding
					n_rho(rr) = n_clad;
				end
			end
			dn_fiber = n_rho - n_clad;
		end

		function [dn_fiber, n_clad, rho_arr] = get_nr_alpha_law_graded_trench_profile(nr, dr, relative_dn, a_core_radius_m, n_clad, trench_width_m, trench_dip, alpha_power)
			% Source: Physical modeling of 10 gbe optical communication systems - Gholami, Molin and Sillard, JLT 2011
			% Trench: Sillard, P., Molin, D., Bigot-Astruc, M., de Jongh, K., & Achten, F. (2017). Rescaled Multimode Fibers for Mode-Division Multiplexing. Journal of Lightwave Technology, 35(8), 1444–1449.
			% dn_fiber = index variation in (x,y) that is an offset from the cladding index, i.e. n(x,y) = dn_fiber(x,y) + n_clad
			rho_arr = linspace(0, nr*dr, nr);

			n_core = n_clad / sqrt(1 - 2*relative_dn);
			n_rho = n_clad * ones(1, nr);
			for rr = 1:nr
				ro = rho_arr(rr);
				if ro <= a_core_radius_m
					% inside the core
					n_rho(rr) = n_core * sqrt(1 - 2*relative_dn*(ro/a_core_radius_m)^alpha_power);
				elseif ro <= a_core_radius_m + trench_width_m
					% inside the trench
					n_rho(rr) = n_clad - trench_dip;
				else
					n_rho(rr) = n_clad;
				end
			end
			dn_fiber = n_rho - n_clad;
		end

		function [dn_fiber, n_clad, x_grid, y_grid] = get_nxy_graded_index_ring_core_profile(nx, ny, dx, dy, relative_dn, ring_radius_m, ring_thickness_m, n_clad, alpha_power)
			% Jin, X. Q., Gomez, A., O’Brien, D. C., & Payne, F. P. (2014). Influence of refractive index profile of ring-core fibres for space division multiplexing systems. Proceedings - 2014 Summer Topicals Meeting Series, SUM 2014, 178–179. https://doi.org/10.1109/SUM.2014.98
			x_grid = linspace(-nx*dx/2, nx*dx/2, nx);
			y_grid = linspace(-ny*dy/2, ny*dy/2, ny);

			n_core = n_clad / sqrt(1 - 2*relative_dn);
			n_xy = n_clad * ones(nx, ny);
			for xx = 1:nx
				for yy = 1:ny
					ro = sqrt(x_grid(xx)^2 + y_grid(yy)^2);
					if abs(ro - ring_radius_m)  <= 0.5*ring_thickness_m
						% inside the ring
						n_xy(yy,xx) = n_core * sqrt(1 - 2*relative_dn*((ro-ring_radius_m)/(ring_thickness_m/2))^alpha_power);
					else
						% outside the ring
						n_xy(yy,xx) = n_clad;
					end
				end
			end
			dn_fiber = n_xy - n_clad;
		end

		function [dn_fiber, n_clad, rho_arr] = get_nr_graded_index_ring_core_profile(nr, dr, relative_dn, ring_radius_m, ring_thickness_m, n_clad, alpha_power)
			% Jin, X. Q., Gomez, A., O’Brien, D. C., & Payne, F. P. (2014). Influence of refractive index profile of ring-core fibres for space division multiplexing systems. Proceedings - 2014 Summer Topicals Meeting Series, SUM 2014, 178–179. https://doi.org/10.1109/SUM.2014.98
			rho_arr = linspace(0, nr*dr, nr);

			n_core = n_clad / sqrt(1 - 2*relative_dn);
			n_rho = n_clad * ones(1, nr);
			for rr = 1:nr
				ro = rho_arr(rr);
				if abs(ro - ring_radius_m)  <= 0.5*ring_thickness_m
					% inside the ring
					n_rho(rr) = n_core * sqrt(1 - 2*relative_dn*((ro-ring_radius_m)/(ring_thickness_m/2))^alpha_power);
				else
					% outside the ring
					n_rho(rr) = n_clad;
				end
			end
			dn_fiber = n_rho - n_clad;
		end


		function [dn_fiber, n_clad, rho_arr] = get_nr_step_index_profile(nr, dr, relative_dn, a_core_radius_m, n_clad)
			% Source: Physical modeling of 10 gbe optical communication systems - Gholami, Molin and Sillard, JLT 2011
			% dn_fiber = index variation in (x,y) that is an offset from the cladding index, i.e. n(x,y) = dn_fiber(x,y) + n_clad
			rho_arr = linspace(0, nr*dr, nr);

			n_core = n_clad / sqrt(1 - 2*relative_dn);
			n_rho = n_clad * ones(1, nr);
			for rr = 1:nr
				ro = rho_arr(rr);
				if ro <= a_core_radius_m
					% inside the core
					n_rho(rr) = n_core;
				else
					% inside the cladding
					n_rho(rr) = n_clad;
				end
			end
			dn_fiber = n_rho - n_clad;
		end



	end
end