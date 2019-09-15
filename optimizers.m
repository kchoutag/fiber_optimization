classdef optimizers
	properties
	end
	methods
		function obj = optimizers()
		end

		function [d_index_distr] = opt_mode_coupling_freeform(obj, fiber_params, init_fiber_params, opt_params)
			D = fiber_params('D');
			nx = fiber_params('nx');
			ny = fiber_params('ny');
			n_clad = utils.get_index_at_wavelength(fiber_params('center_wavelength_nm'));

			init_neff = init_fiber_params('neff');
			curr_neff = fiber_params('neff');
			curr_fields = fiber_params('fields');

			% set the optimization target
			if(opt_params('direction') == 'MIN')
				des_neff = linspace(init_neff(1), init_neff(end), init_fiber_params('D'));
			elseif (opt_params('direction') == 'MAX')
				des_neff = ones(1,init_fiber_params('D'))*mean(init_neff);
            end

            % compute the index update
			d_index_distr = zeros(ny,nx);
			for xx = 1:nx
				for yy = 1:ny
					for mm = 1:min(init_fiber_params('D'), D)
						d_index_distr(yy,xx) = d_index_distr(yy,xx) - ((curr_neff(mm) - des_neff(mm))/(curr_neff(mm)))*curr_fields(yy,xx,mm)^2;
					end
				end
			end
			d_index_distr = d_index_distr .* (fiber_params('nxy_offset_from_cladding') + n_clad);

			% normalize index update
			d_index_distr = opt_params('max_dn')*d_index_distr/max(max(d_index_distr));
		end

		function [d_n_rho] = opt_mode_coupling_radial(obj, fiber_params, init_fiber_params, opt_params)
			D = fiber_params('D');
			nr = fiber_params('nr');
			dr = fiber_params('dr');
			n_clad = utils.get_index_at_wavelength(fiber_params('center_wavelength_nm'));

			init_neff = init_fiber_params('neff');
			curr_neff = fiber_params('neff');
			curr_fields = fiber_params('fields');

			rho_arr = linspace(0, nr*dr, nr);
			x_arr = linspace(-nr*dr, nr*dr, nr);
			y_arr = x_arr;

			% set the optimization target
			if(opt_params('direction') == 'MIN')
				des_neff = linspace(init_neff(1), init_neff(end), init_fiber_params('D'));
			elseif (opt_params('direction') == 'MAX')
				des_neff = ones(1,init_fiber_params('D'))*mean(init_neff);
            end

            % compute the index update
            d_n_rho = zeros(1, nr);
            for rr = 1:nr
            	rho = rho_arr(rr);
            	d_n_rho(rr) = 0;
            	for mm = 1: min(init_fiber_params('D'), D)
            		mode_mm_intensity_at_rho = utils.integrate_mode_on_circle(curr_fields(:,:,mm), x_arr, y_arr, rho);
            		d_n_rho(rr) = d_n_rho(rr) - ((curr_neff(mm) - des_neff(mm))/curr_neff(mm))*mode_mm_intensity_at_rho;
            	end
            end
            d_n_rho = d_n_rho .* (fiber_params('nr_offset_from_cladding') + n_clad);

            % normalize index update
            d_n_rho = opt_params('max_dn')*d_n_rho/max(d_n_rho);
		end


	end
end