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
			opt_params('opt_steps') = 100;
			opt_params('max_dn') = 1e-4;
			opt_params('direction') = "MIN";

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
				u.plot_cell_array('neff evolution', neff_hist, 'Iteration', 'Effective Index');
				u.plot_results('Fiber', fiber_params);

				drawnow;
			end
		end

	end
end