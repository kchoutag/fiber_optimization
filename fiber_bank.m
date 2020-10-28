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

		function [fib_params] = get_MJL_Alpha2p095(obj, center_wavelength_nm, B) % graded index alpha=2.095 profile (axially-symmetric)
			% B is a list of coefficients from MJL's patent

			u = utils();
            fib_params = containers.Map;

			fib_params('nr') = 100;
			fib_params('dr') = 0.3e-6;
			fib_params('center_wavelength_nm') = center_wavelength_nm;
			fib_params('d_wavelength_nm') = 0.001;
			fib_params('relative_index_diff_delta') = 0.01;
			fib_params('a_core_radius_m') = 25*1e-6;
			fib_params('alpha_power') = 2.095; 
			fib_params('axially_symm') = true;
			fib_params('D_upper_limit') = 120;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, rho_arr, alpha_part, non_alpha_part] = obj.profile_generator.get_nr_MJL_profile(B, fib_params('nr'), fib_params('dr'), fib_params('relative_index_diff_delta'), fib_params('a_core_radius_m'), n_clad, fib_params('alpha_power'));
			
			fib_params('nr_offset_from_cladding') = dn_fiber;

			fib_params('nr_alpha_part') 	= alpha_part;
			fib_params('nr_non_alpha_part') = non_alpha_part;
		end

		function [fib_params] = get_GI_MMF_2(obj) % graded index alpha=2 profile (axially-symmetric)
			u = utils();
            fib_params = containers.Map;

			fib_params('nr') = 100;
			fib_params('dr') = 0.25e-6;
			fib_params('center_wavelength_nm') = 1550;
			fib_params('d_wavelength_nm') = 0.001;
			fib_params('relative_index_diff_delta') = 0.01;
			fib_params('a_core_radius_m') = 15*1e-6;
			fib_params('alpha_power') = 2;
			fib_params('axially_symm') = true;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~] = obj.profile_generator.get_nr_alpha_law_graded_profile(fib_params('nr'), fib_params('dr'), fib_params('relative_index_diff_delta'), fib_params('a_core_radius_m'), n_clad, fib_params('alpha_power'));
			fib_params('nr_offset_from_cladding') = dn_fiber;
		end

		function [fib_params] = get_StepIndex_MMF(obj) % (axially-symmetric)
			% Fiber from: P Sillard et al. Few-mode fiber for uncoupled mode-division multiplexing transmissions. ECOC 2011, 38–40.
			u = utils();
            fib_params = containers.Map;

			fib_params('nr') = 100;
			fib_params('dr') = 0.20e-6;
			fib_params('center_wavelength_nm') = 1550; 
			fib_params('d_wavelength_nm') = 0.001;
			fib_params('relative_index_diff_delta') = 0.00665; % need n_co - n_cl = 9.7*1e-3 according to the paper
			fib_params('a_core_radius_m') = 7.5*1e-6;
			fib_params('axially_symm') = true;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~] = obj.profile_generator.get_nr_step_index_profile(fib_params('nr'), fib_params('dr'), fib_params('relative_index_diff_delta'), fib_params('a_core_radius_m'), n_clad);
			fib_params('nr_offset_from_cladding') = dn_fiber;
		end

		function [fib_params] = get_StepIndex_MMF_debug(obj) % (axially-symmetric)
			u = utils();
            fib_params = containers.Map;

			fib_params('nr') = 100;
			fib_params('dr') = 0.40e-6;
			fib_params('center_wavelength_nm') = 1550; 
			fib_params('d_wavelength_nm') = 0.001;
			fib_params('relative_index_diff_delta') = 0.00665;
			fib_params('a_core_radius_m') = 20*1e-6;
			fib_params('axially_symm') = true;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~] = obj.profile_generator.get_nr_step_index_profile(fib_params('nr'), fib_params('dr'), fib_params('relative_index_diff_delta'), fib_params('a_core_radius_m'), n_clad);
			fib_params('nr_offset_from_cladding') = dn_fiber;
		end

		function [fib_params] = get_RCF_1(obj) % graded-index ring-core fiber
			% Feng, F., Guo, X., Gordon, G. S. D., Jin, X. Q., Payne, F. P., Jung, Y., … Wilkinson, T. D. (2016). All-optical mode-group division multiplexing over a graded-index ring-core fiber with single radial mode. 2016 Optical Fiber Communications Conference and Exhibition, OFC 2016, 3–5.
			u = utils();
            fib_params = containers.Map;

			fib_params('nr') = 100;
			fib_params('dr') = 0.20e-6;
			fib_params('center_wavelength_nm') = 1550;
			fib_params('d_wavelength_nm') = 0.001;
			fib_params('relative_index_diff_delta') = 0.01;
			fib_params('ring_radius_m') = 10*1e-6;
			fib_params('ring_thickness_m') = 4.3*1e-6;
			fib_params('alpha_power') = 8;
			fib_params('axially_symm') = true;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~] = obj.profile_generator.get_nr_graded_index_ring_core_profile(fib_params('nr'), fib_params('dr'), fib_params('relative_index_diff_delta'), fib_params('ring_radius_m'), fib_params('ring_thickness_m'), n_clad, fib_params('alpha_power'));
			fib_params('nr_offset_from_cladding') = dn_fiber;
		end

		function [fib_params] = get_RCF_2(obj) % graded-index ring-core fiber (free-form)
			% Parameters different from: Feng, F., Guo, X., Gordon, G. S. D., Jin, X. Q., Payne, F. P., Jung, Y., … Wilkinson, T. D. (2016). All-optical mode-group division multiplexing over a graded-index ring-core fiber with single radial mode. 2016 Optical Fiber Communications Conference and Exhibition, OFC 2016, 3–5.
			u = utils();
            fib_params = containers.Map;

			fib_params('nx') = 100;
			fib_params('ny') = 100;
			fib_params('dx') = 0.20e-6;
			fib_params('dy') = 0.20e-6;
			fib_params('center_wavelength_nm') = 1550;
			fib_params('d_wavelength_nm') = 0.001;
			fib_params('relative_index_diff_delta') = 0.01;
			fib_params('ring_radius_m') = 6*1e-6;
			fib_params('ring_thickness_m') = 4*1e-6;
			fib_params('alpha_power') = 8;
			fib_params('axially_symm') = false;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~] = obj.profile_generator.get_nxy_graded_index_ring_core_profile(fib_params('nx'), fib_params('ny'), fib_params('dx'), fib_params('dy'), fib_params('relative_index_diff_delta'), fib_params('ring_radius_m'), fib_params('ring_thickness_m'), n_clad, fib_params('alpha_power'));
			fib_params('nxy_offset_from_cladding') = dn_fiber;
		end

		function [fib_params] = get_RCF_3(obj) % ring-core fiber (free-form)
			% Parameters from: Kasahara, M., Saitoh, K., Sakamoto, T., Hanzawa, N., Matsui, T., Tsujikawa, K., & Yamamoto, F. (2014). Design of three-spatial-mode ring-core fiber. Journal of Lightwave Technology, 32(7), 1337–1343. https://doi.org/10.1109/JLT.2014.2304732
			u = utils();
            fib_params = containers.Map;

			fib_params('nx') = 100;
			fib_params('ny') = 100;
			fib_params('dx') = 0.20e-6;
			fib_params('dy') = 0.20e-6;
			fib_params('center_wavelength_nm') = 1550;
			fib_params('d_wavelength_nm') = 0.001;
			fib_params('relative_index_diff_delta') = 0.8*1e-2;
			a = 7*1e-6; d = 2*a*0.3; % parameters from Kasahara, M et al. for d/2a = 0.3
			fib_params('ring_radius_m') = a/2 + d/4; % converting to notation from X Jin et al.
			fib_params('ring_thickness_m') = a - d/2; % converting to X Jin et al.
			fib_params('alpha_power') = 100; % high alpha = effectively step index
			fib_params('axially_symm') = false;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~] = obj.profile_generator.get_nxy_graded_index_ring_core_profile(fib_params('nx'), fib_params('ny'), fib_params('dx'), fib_params('dy'), fib_params('relative_index_diff_delta'), fib_params('ring_radius_m'), fib_params('ring_thickness_m'), n_clad, fib_params('alpha_power'));
			fib_params('nxy_offset_from_cladding') = dn_fiber;
		end

		function [fib_params] = get_4C_MCF_debug(obj)
			u = utils();
            fib_params = containers.Map;

            fib_params('nx') = 100;
			fib_params('ny') = 100;
			fib_params('dx') = 0.30e-6;
			fib_params('dy') = 0.30e-6;
			fib_params('center_wavelength_nm') = 1550;
			fib_params('d_wavelength_nm') = 0.01;
			fib_params('relative_index_diff_delta') = 0.01;
			fib_params('n_cores') = 4;
			fib_params('core_radii_m') = 3*1e-6 *ones(1,fib_params('n_cores'));
			fib_params('core_xlocs_m') = [1 1 -1 -1]*4*1e-6;
			fib_params('core_ylocs_m') = [1 -1 1 -1]*4*1e-6;

			fib_params('axially_symm') = false;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~, ~] = obj.profile_generator.get_nxy_multicore_profile(fib_params('nx'), fib_params('ny'), fib_params('dx'), fib_params('dy'), fib_params('n_cores'), fib_params('core_radii_m'), fib_params('core_xlocs_m'), fib_params('core_ylocs_m'), n_clad, fib_params('relative_index_diff_delta'));
			fib_params('nxy_offset_from_cladding') = dn_fiber;
		end


		function [fib_params] = get_MMF_withTrench_1(obj) % 15-mode strongly-coupled MMF from Sillard et al JLT 2016 Low-Differential-Mode-Group-Delay 9-LP-Mode Fiber
			u = utils();
            fib_params = containers.Map;

			fib_params('nr') = 100;
			fib_params('dr') = 0.25e-6;
			fib_params('center_wavelength_nm') = 1550;
			fib_params('a_core_radius_m') = 14*1e-6;
			fib_params('d_wavelength_nm') = 0.01;
			fib_params('relative_index_diff_delta') = 0.0102; % chosen to get a max core-cladding difference of 14.9*1e-3 as per paper.
			fib_params('trench_width_m') = 4*1e-6;
			fib_params('trench_dip') = 6*1e-3;
			fib_params('alpha_power') = 1.94;
			fib_params('axially_symm') = true;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~] = obj.profile_generator.get_nr_alpha_law_graded_trench_profile(fib_params('nr'), fib_params('dr'), fib_params('relative_index_diff_delta'), fib_params('a_core_radius_m'), n_clad, fib_params('trench_width_m'), fib_params('trench_dip'), fib_params('alpha_power'));
			fib_params('nr_offset_from_cladding') = dn_fiber;
		end

		function [fib_params] = get_Corning_SMF28(obj) % step index profile (axially-symmetric)
			%Fiber from Riishede, J., & Sigmund, O. (2008). Inverse design of dispersion compensating optical fiber using topology optimization. Journal of the Optical Society of America B: Optical Physics, 25(1), 88–97. https://doi.org/10.1364/JOSAB.25.000088
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

		function [fib_params] = get_SMF_debug(obj)
			u = utils();
            fib_params = containers.Map;

			fib_params('nr') = 100;
			fib_params('dr') = 0.20e-6;
			fib_params('center_wavelength_nm') = 1550;
			fib_params('d_wavelength_nm') = 0.01;
			fib_params('relative_index_diff_delta') = .2 * 1e-2;
			fib_params('a_core_radius_m') = 7*1e-6;
			fib_params('axially_symm') = true;

			n_clad = u.get_index_at_wavelength(fib_params('center_wavelength_nm'));
			[dn_fiber, ~, ~] = obj.profile_generator.get_nr_step_index_profile(fib_params('nr'), fib_params('dr'), fib_params('relative_index_diff_delta'), fib_params('a_core_radius_m'), n_clad);
			fib_params('nr_offset_from_cladding') = dn_fiber;
		end

	end
end