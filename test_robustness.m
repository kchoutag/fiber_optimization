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

			rms_gd_variation = {};
			for ii = 1:length(lambda_arr)
				disp(sprintf('Checking wavelength #%d of %d...', ii, length(lambda_arr)));
				copy_fib_params('center_wavelength_nm') = lambda_arr(ii)*1e9;
				copy_fib_params = u.solve_fiber_properties(copy_fib_params);
				rms_gd_variation{end + 1} = rms(copy_fib_params('MD_coeffs_psm'));
			end
		end

	end
end