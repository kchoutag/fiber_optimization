classdef fiber_bank
	% Bank of pre-designed fibers
	properties
		profile_generator = profiles();
		c_speed_light = 299792458; % meters per second
	end

	methods
		function obj = fiber_bank();
		end

		function [fib_params] = get_GI_MMF_1(obj) % graded index / alpha law profile
			fib_params = containers.Map;

			fib_params('nx') = 100;
			fib_params('ny') = 100;
			fib_params('dx') = 0.75e-6;
			fib_params('dy') = 0.75e-6;
			fib_params('center_wavelength_nm') = 1550;
			fib_params('d_wavelength_nm') = 0.001;
			fib_params('relative_index_diff_delta') = 0.01;
			fib_params('a_core_radius_m') = 20*1e-6;
			fib_params('n0_core') = 1.44;
			fib_params('alpha_power') = 2;
			fib_params('index_distr') = obj.profile_generator.get_alpha_law_graded_profile(fib_params('nx'), fib_params('ny'), fib_params('dx'), fib_params('dy'), fib_params('relative_index_diff_delta'), fib_params('a_core_radius_m'), fib_params('n0_core'), fib_params('alpha_power'));
		end


		function [fib_params] = solve_fiber_properties(obj, fib_params)
			addpath('../modesolver-2011-04-22/');
			n_modes_upper_lim = 70;
			eps_profile = fib_params('index_distr').^2;
			n_core = max(max(fib_params('index_distr'))); 
			n_cladding = min(min(fib_params('index_distr'))); 

			% center wavelength
			[fields, neff_aux] = svmodes(fib_params('center_wavelength_nm')*1e-9, n_core, n_modes_upper_lim, fib_params('dx'), fib_params('dy'), eps_profile, '0000', 'scalar');
			valid_idx = neff_aux > n_cladding;
			fields = fields(:,:,valid_idx);
			neff = neff_aux(valid_idx);
			D_center = length(neff);

			% left wavelength
			[fields_left, neff_aux_left] = svmodes((fib_params('center_wavelength_nm') - fib_params('d_wavelength_nm'))*1e-9, n_core, n_modes_upper_lim, fib_params('dx'), fib_params('dy'), eps_profile, '0000', 'scalar');
			valid_idx = neff_aux_left > n_cladding;
			fields_left = fields_left(:,:,valid_idx);
			neff_left = neff_aux_left(valid_idx);
			D_left = length(neff_left);

			% right wavelength
			[fields_right, neff_aux_right] = svmodes((fib_params('center_wavelength_nm') + fib_params('d_wavelength_nm'))*1e-9, n_core, n_modes_upper_lim, fib_params('dx'), fib_params('dy'), eps_profile, '0000', 'scalar');
			valid_idx = neff_aux_right > n_cladding;
			fields_right = fields_right(:,:,valid_idx);
			neff_right = neff_aux_right(valid_idx);
			D_right = length(neff_right);

			if ((D_center ~= D_left) || (D_center ~= D_right) || (D_left ~= D_right))
				error('Number of modes not equal at left, center, right wavelengths!');
			end

			D = min(min(D_center, D_left), D_right);
			
			% calculate neff = beta0 of modes = phase velocities
			neff = neff(1:D);
			neff_left = neff_left(1:D);
			neff_right = neff_right(1:D);

			% set the fields
			fields = fields(:,:,1:D);
			fields_left = fields_left(:,:,1:D);
			fields_right = fields_right(:,:,1:D);

			% calculate the modal dispersion for all modes using the central difference approximation (in arbitrary units, TODO: fix units?)
    		modal_dispersion_coeffs = (neff_right - neff_left)/(2*fib_params('d_wavelength_nm')*1e-9);

    		% calculate chromatic dispersion for all modes (in arbitrary units)
    		chromatic_dispersion_coeffs = (-fib_params('center_wavelength_nm')*1e-9/(obj.c_speed_light * (fib_params('d_wavelength_nm')*1e-9)^2)) * ...
			    				  (neff_left + neff_right - 2*neff);

			% store all parameters
			fib_params('D') = D;
			fib_params('fields') = fields; fib_params('fields_left') = fields_left; fib_params('fields_right') = fields_right;
			fib_params('neff') = neff; fib_params('neff_left') = neff_left; fib_params('neff_right') = neff_right;

			fib_params('modal_dispersion_coeffs') = modal_dispersion_coeffs;
			fib_params('chromatic_dispersion_coeffs') = chromatic_dispersion_coeffs;
			
		end

	end
end