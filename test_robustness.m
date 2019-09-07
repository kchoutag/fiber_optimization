classdef test_robustness
	methods
		function obj = test_robustness
		end

		function [neff_variation] = vary_wavelength_neff(obj, fiber_params, lambda_arr)
			b = fiber_bank();
			copy_fib_params = containers.Map(fiber_params.keys, fiber_params.values);

			neff_variation = {};
			for ii = 1:length(lambda_arr)
				disp(sprintf('Checking wavelength #%d of %d...', ii, length(lambda_arr)));
				copy_fib_params('center_wavelength_nm') = lambda_arr(ii)*1e9;
				copy_fib_params = b.solve_fiber_properties(copy_fib_params);
				neff_variation{end + 1} = copy_fib_params('neff');
			end

		end

	end
end