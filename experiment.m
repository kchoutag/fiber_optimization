classdef experiment
	properties
		bank = fiber_bank();
		
		optimizer = optimizers();
	end
	methods
		function obj = experiment()
		end

		function obj = run_MMF_GI_reduce_coupling(obj)
			fiber_params = obj.bank.get_GI_MMF_1();
			fiber_params = obj.bank.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 80;
			opt_params('max_dn') = 1e-4;
			opt_params('direction') =  "MIN"; %"MAX"; %"MIN";

			neff_hist = {fiber_params('neff')};

			u = utils();

			for nn = 1:opt_params('opt_steps')
				disp(sprintf('Running iteration %d of %d...', nn, opt_params('opt_steps')));
				% compute and apply the index update
				index_update = obj.optimizer.opt_mode_coupling_freeform(fiber_params, init_fiber_params, opt_params);
				fiber_params('index_distr') = fiber_params('index_distr') + index_update;

				% update the realized fiber properties
				fiber_params = obj.bank.solve_fiber_properties(fiber_params);

				% update the history
				neff_hist{end+1} = fiber_params('neff');
				u.plot_cell_array('neff evolution', 1:nn+1, neff_hist, 'Iteration', 'Effective Index');
				u.plot_results('Fiber', fiber_params);

				drawnow;
			end

			robustness = test_robustness();
			lambda_arr = linspace(1200,1700,20)*1e-9;
			neff_wavelength_variation = robustness.vary_wavelength_neff(fiber_params, lambda_arr);
			u.plot_cell_array('wavelength variation', lambda_arr*1e9, neff_wavelength_variation, 'Wavelength (nm)', 'Effective Index');
			axis tight;
		end

		function obj = run_playground(obj) % workspace for debugging and/or developing new fiber designs
			fiber_params = obj.bank.get_GI_MMF_1();
			fiber_params = obj.bank.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 2;
			opt_params('max_dn') = 1e-4;
			opt_params('direction') =  "MIN"; %"MAX"; %"MIN";

			neff_hist = {fiber_params('neff')};

			u = utils();

			for nn = 1:opt_params('opt_steps')
				disp(sprintf('Running iteration %d of %d...', nn, opt_params('opt_steps')));
				% compute and apply the index update
				index_update = obj.optimizer.opt_mode_coupling_freeform(fiber_params, init_fiber_params, opt_params);
				fiber_params('index_distr') = fiber_params('index_distr') + index_update;

				% update the realized fiber properties
				fiber_params = obj.bank.solve_fiber_properties(fiber_params);

				% update the history
				neff_hist{end+1} = fiber_params('neff');
				u.plot_cell_array('neff evolution', 1:nn+1, neff_hist, 'Iteration', 'Effective Index');
				u.plot_results('Fiber', fiber_params);

				drawnow;
			end

			robustness = test_robustness();
			lambda_arr = linspace(1200,1700,20)*1e-9;
			neff_wavelength_variation = robustness.vary_wavelength_neff(fiber_params, lambda_arr);
			u.plot_cell_array('wavelength variation', lambda_arr*1e9, neff_wavelength_variation, 'Wavelength (nm)', 'Effective Index');
			axis tight;
		end

	end
end