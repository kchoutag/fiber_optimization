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

			D_optimize = min(init_fiber_params('D'), D);

			% set the optimization target
			if(strcmp(opt_params('direction'), 'MIN'))
				des_neff = linspace(init_neff(1), init_neff(D_optimize), D_optimize);
			elseif(strcmp(opt_params('direction'),'MAX'))
				des_neff = ones(1,D_optimize)*mean(init_neff);
			elseif(strcmp(opt_params('direction'),'RCF_to_MCF'))
				%Target the same core index but reduce the current spread of neff during each iteration
				des_neff = init_neff(1) + 0.5*linspace(0, curr_neff(D_optimize)-curr_neff(1), D_optimize);
            end

            % compute the index update
			d_index_distr = zeros(ny,nx);
			for xx = 1:nx
				for yy = 1:ny
					for mm = 1:D_optimize
						d_index_distr(yy,xx) = d_index_distr(yy,xx) - ((curr_neff(mm) - des_neff(mm))/(curr_neff(mm)))*curr_fields(yy,xx,mm)^2;
					end
				end
			end
			d_index_distr = d_index_distr .* (fiber_params('nxy_offset_from_cladding') + n_clad);

			% normalize index update
			d_index_distr = opt_params('max_dn')*d_index_distr/max(max(abs(d_index_distr)));
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
			if(strcmp(opt_params('direction'), 'MIN'))
				des_neff = linspace(init_neff(1), init_neff(end), init_fiber_params('D'));
			elseif(strcmp(opt_params('direction'),'MAX'))
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
            small_val = 1e-12; % to prevent divide by zero
            d_n_rho = opt_params('max_dn')*d_n_rho/max(abs(d_n_rho) + small_val);
		end

		function [d_n_rho] = opt_modal_dispersion_radial(obj, fiber_params, init_fiber_params, opt_params)
			D = fiber_params('D');
			nr = fiber_params('nr');
			dr = fiber_params('dr');
			rho_arr = linspace(0, nr*dr, nr);
			x_arr = linspace(-nr*dr, nr*dr, nr);
			y_arr = x_arr;

			n_clad = utils.get_index_at_wavelength(fiber_params('center_wavelength_nm'));

			neff_left = fiber_params('neff_left'); neff_right = fiber_params('neff_right');
			fields_left = fiber_params('fields_left'); fields_right = fiber_params('fields_right'); 

			MD_coeffs_psm = fiber_params('MD_coeffs_psm'); % MD coefficients (= group delay ps per unit m) are already zero-mean

			% set the optimization target
			if(strcmp(opt_params('direction'),'MIN')) %TODO: check the signs??
				step_sign = 1;
			else
				step_sign = -1;
            end

            % compute the index update
            n_rho_curr = (fiber_params('nr_offset_from_cladding') + n_clad);
            d_n_rho = zeros(1, nr);
            for rr = 1:nr
            	rho = rho_arr(rr);
            	d_n_rho(rr) = 0;
            	for mm = 1: min(init_fiber_params('D'), D)
            		mode_mm_left_intensity_at_rho   = utils.integrate_mode_on_circle(fields_left(:,:,mm), x_arr, y_arr, rho);
            		mode_mm_right_intensity_at_rho  = utils.integrate_mode_on_circle(fields_right(:,:,mm), x_arr, y_arr, rho);

            		d_MDi_dn = -1*n_rho_curr(rr)*( mode_mm_right_intensity_at_rho/neff_right(mm)...
            									 - mode_mm_left_intensity_at_rho/neff_left(mm));

            		% minimize the rms MD --> d(rms MD)/dn proportional to sum_i(MD_i * dMD_i/dn)
            		d_n_rho(rr) = d_n_rho(rr) + MD_coeffs_psm(mm)*d_MDi_dn;
            	end
            end

            % normalize index update
            d_n_rho = step_sign*opt_params('max_dn')*d_n_rho/max(abs(d_n_rho));
		end

		function [d_n_rho] = opt_chromatic_dispersion_radial(obj, fiber_params, init_fiber_params, opt_params)
			D = fiber_params('D');
			nr = fiber_params('nr');
			dr = fiber_params('dr');
			rho_arr = linspace(0, nr*dr, nr);
			x_arr = linspace(-nr*dr, nr*dr, nr);
			y_arr = x_arr;

			n_clad = utils.get_index_at_wavelength(fiber_params('center_wavelength_nm'));

			neff_center = fiber_params('neff'); neff_left = fiber_params('neff_left'); neff_right = fiber_params('neff_right');
			fields_center = fiber_params('fields'); fields_left = fiber_params('fields_left'); fields_right = fiber_params('fields_right'); 

			CD_coeffs_psnmkm = fiber_params('CD_coeffs_psnmkm');

			% set the optimization target
			if(strcmp(opt_params('direction'),'MIN'))
				step_sign = -1;
			else
				step_sign = 1;
            end

            % compute the index update
            n_rho_curr = (fiber_params('nr_offset_from_cladding') + n_clad);
            d_n_rho = zeros(1, nr);
            for rr = 1:nr
            	rho = rho_arr(rr);
            	d_n_rho(rr) = 0;
            	for mm = 1: min(init_fiber_params('D'), D)
            		mode_mm_center_intensity_at_rho = utils.integrate_mode_on_circle(fields_center(:,:,mm), x_arr, y_arr, rho);
            		mode_mm_left_intensity_at_rho   = utils.integrate_mode_on_circle(fields_left(:,:,mm), x_arr, y_arr, rho);
            		mode_mm_right_intensity_at_rho  = utils.integrate_mode_on_circle(fields_right(:,:,mm), x_arr, y_arr, rho);

            		d_CDi_dn = -1*n_rho_curr(rr)*( mode_mm_left_intensity_at_rho/neff_left(mm)  ...
						            			 + mode_mm_right_intensity_at_rho/neff_right(mm) ...
						            			 - 2*mode_mm_center_intensity_at_rho/neff_center(mm));

            		% minimize the rms CD --> d(rms CD)/dn proportional to sum_i(CD_i * dCD_i/dn)
            		d_n_rho(rr) = d_n_rho(rr) + CD_coeffs_psnmkm(mm)*d_CDi_dn;
            	end
            end

            % normalize index update
            d_n_rho = step_sign*opt_params('max_dn')*d_n_rho/max(abs(d_n_rho));
		end

	end
end