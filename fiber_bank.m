classdef fiber_bank
	% Bank of pre-designed fibers
	properties
		profile_generator = profiles();
	end

	methods
		function obj = fiber_bank()
		end

		function [fib_params] = get_GI_MMF_1(obj) % graded index alpha=2 profile (free-form)
			u = utils();
            fib_params = containers.Map;

			fib_params('nx') = 100;
			fib_params('ny') = 100;
			fib_params('dx') = 0.3e-6;
			fib_params('dy') = 0.3e-6;
			fib_params('center_wavelength_nm') = 1550;
			fib_params('d_wavelength_nm') = 0.01;
			fib_params('relative_index_diff_delta') = 0.01;
			fib_params('a_core_radius_m') = 10*1e-6;
			fib_params('alpha_power') = 2;
			fib_params('axially_symm') = false;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~, ~] = obj.profile_generator.get_nxy_alpha_law_graded_profile(fib_params('nx'), fib_params('ny'), fib_params('dx'), fib_params('dy'), fib_params('relative_index_diff_delta'), fib_params('a_core_radius_m'), n_clad, fib_params('alpha_power'));
			fib_params('nxy_offset_from_cladding') = dn_fiber;
		end

		function [fib_params] = get_GI_MMF_2(obj) % graded index alpha=2 profile (axially-symmetric)
			u = utils();
            fib_params = containers.Map;

			fib_params('nr') = 100;
			fib_params('dr') = 0.15e-6;
			fib_params('center_wavelength_nm') = 1550;
			fib_params('d_wavelength_nm') = 0.01;
			fib_params('relative_index_diff_delta') = 0.01;
			fib_params('a_core_radius_m') = 10*1e-6;
			fib_params('alpha_power') = 2;
			fib_params('axially_symm') = true;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~] = obj.profile_generator.get_nr_alpha_law_graded_profile(fib_params('nr'), fib_params('dr'), fib_params('relative_index_diff_delta'), fib_params('a_core_radius_m'), n_clad, fib_params('alpha_power'));
			fib_params('nr_offset_from_cladding') = dn_fiber;
		end

		function [fib_params] = get_GI_SMF_1(obj) % graded index alpha=2 profile (axially-symmetric)
			u = utils();
            fib_params = containers.Map;

			fib_params('nr') = 100;
			fib_params('dr') = 0.07e-6;
			fib_params('center_wavelength_nm') = 1550;
			fib_params('d_wavelength_nm') = 0.01;
			fib_params('relative_index_diff_delta') = 0.01;
			fib_params('a_core_radius_m') = 2*1e-6;
			fib_params('alpha_power') = 2;
			fib_params('axially_symm') = true;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~] = obj.profile_generator.get_nr_alpha_law_graded_profile(fib_params('nr'), fib_params('dr'), fib_params('relative_index_diff_delta'), fib_params('a_core_radius_m'), n_clad, fib_params('alpha_power'));
			fib_params('nr_offset_from_cladding') = dn_fiber;
		end

		function [fib_params] = get_SI_SMF_1(obj) % step index profile (axially-symmetric)
			u = utils();
            fib_params = containers.Map;

			fib_params('nr') = 100;
			fib_params('dr') = 0.50e-6;
			fib_params('center_wavelength_nm') = 1550;
			fib_params('d_wavelength_nm') = 0.01;
			fib_params('relative_index_diff_delta') = .36 * 1e-2;
			fib_params('a_core_radius_m') = 4.1*1e-6;
			fib_params('axially_symm') = true;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~] = obj.profile_generator.get_nr_step_index_profile(fib_params('nr'), fib_params('dr'), fib_params('relative_index_diff_delta'), fib_params('a_core_radius_m'), n_clad);
			fib_params('nr_offset_from_cladding') = dn_fiber;
		end

	end
end