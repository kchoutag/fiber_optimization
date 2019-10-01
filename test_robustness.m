classdef test_robustness
	methods(Static)
		function [neff_variation] = vary_wavelength_neff(fiber_params, lambda_arr)
            u = utils();
			copy_fib_params = containers.Map(fiber_params.keys, fiber_params.values);

			neff_variation = {};
			for ii = 1:length(lambda_arr)
				disp(sprintf('Checking wavelength #%d of %d...', ii, length(lambda_arr)));
				copy_fib_params('center_wavelength_nm') = lambda_arr(ii)*1e9;
				copy_fib_params = u.solve_fiber_properties(copy_fib_params);
				neff_variation{end + 1} = copy_fib_params('neff');
			end
		end

		function [rms_gd_variation] = vary_wavelength_rms_gd(fiber_params, lambda_arr)
            u = utils();
			copy_fib_params = containers.Map(fiber_params.keys, fiber_params.values);

			rms_gd_variation = [];
			for ii = 1:length(lambda_arr)
				disp(sprintf('Checking wavelength #%d of %d...', ii, length(lambda_arr)));
				copy_fib_params('center_wavelength_nm') = lambda_arr(ii)*1e9;
				copy_fib_params = u.solve_fiber_properties(copy_fib_params);
				rms_gd_variation = [rms_gd_variation rms(copy_fib_params('MD_coeffs_psm'))];
			end
		end

		function [D_variation] = vary_wavelength_number_modes(fiber_params, lambda_arr)
            u = utils();
			copy_fib_params = containers.Map(fiber_params.keys, fiber_params.values);

			D_variation = [];
			for ii = 1:length(lambda_arr)
				disp(sprintf('Checking wavelength #%d of %d...', ii, length(lambda_arr)));
				copy_fib_params('center_wavelength_nm') = lambda_arr(ii)*1e9;
				copy_fib_params = u.solve_fiber_properties(copy_fib_params);
				D_variation = [D_variation copy_fib_params('D')];
			end
		end

		function [fiber_params_out] = change_gridding(fiber_params, interp_factor)
			u = utils();
			fiber_params_out = containers.Map(fiber_params.keys, fiber_params.values);
			rho_in = linspace(0, fiber_params_out('nr')*fiber_params_out('dr'), fiber_params_out('nr'));
			[fiber_params_out('nr_offset_from_cladding'), rho_out] = u.interpolate_index_radial(fiber_params_out('nr_offset_from_cladding'), rho_in, interp_factor);
			fiber_params_out('nr') = length(fiber_params_out('nr_offset_from_cladding'));
			fiber_params_out('dr') = rho_out(2) - rho_in(1);
		end

	end
end