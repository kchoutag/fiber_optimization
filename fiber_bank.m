classdef fiber_bank
	% Bank of pre-designed fibers
	properties
		profile_generator = profiles();
	end

	methods
		function obj = fiber_bank()
		end

		function [fib_params] = get_GI_MMF_1(obj) % graded index / alpha law profile
			u = utils();
            fib_params = containers.Map;

			fib_params('nx') = 100;
			fib_params('ny') = 100;
			fib_params('dx') = 0.3e-6;
			fib_params('dy') = 0.3e-6;
			fib_params('center_wavelength_nm') = 1550;
			fib_params('d_wavelength_nm') = 0.001;
			fib_params('relative_index_diff_delta') = 0.01;
			fib_params('a_core_radius_m') = 10*1e-6;
			fib_params('alpha_power') = 2;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~, ~] = obj.profile_generator.get_alpha_law_graded_profile(fib_params('nx'), fib_params('ny'), fib_params('dx'), fib_params('dy'), fib_params('relative_index_diff_delta'), fib_params('a_core_radius_m'), n_clad, fib_params('alpha_power'));
			fib_params('index_distr_offset_from_cladding') = dn_fiber;
		end

	end
end