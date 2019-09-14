classdef profiles
	methods(Static)
		function [dn_fiber, n_clad, x_grid, y_grid] = get_alpha_law_graded_profile(nx, ny, dx, dy, relative_dn, a_core_radius_m, n_clad, alpha_power)
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
	end
end