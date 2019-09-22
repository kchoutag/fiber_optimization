classdef experiment
	properties
		bank = fiber_bank();
		
		optimizer = optimizers();
	end
	methods
		function obj = experiment()
		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%           WORKING         %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		function obj = MMF_GI_radial_reduce_modal_dispersion(obj)
			fiber_params = obj.bank.get_MMF_withTrench_1();
			fiber_params = utils.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 100;%80;
			opt_params('max_dn') = 1e-5;
			opt_params('direction') =  "MIN"; 

			rms_MD_hist = [rms(fiber_params('MD_coeffs_psm'))];
			MD_coeffs_hist = {fiber_params('MD_coeffs_psm')};

			for nn = 1:opt_params('opt_steps')
				disp(sprintf('Running iteration %d of %d...', nn, opt_params('opt_steps')));
				% compute and apply the index update
				index_update = obj.optimizer.opt_modal_dispersion_radial(fiber_params, init_fiber_params, opt_params);
				fiber_params('nr_offset_from_cladding') = fiber_params('nr_offset_from_cladding') + index_update;

				% update the realized fiber properties
				fiber_params = utils.solve_fiber_properties(fiber_params);

				% update the history
				rms_MD_hist = [rms_MD_hist rms(fiber_params('MD_coeffs_psm'))];
				MD_coeffs_hist{end+1} = fiber_params('MD_coeffs_psm');

				
				% plot the results
				dsfig('JLT figure 1');
				rho_arr = linspace(0, fiber_params('nr')*fiber_params('dr'), fiber_params('nr'));
	        	n_clad = utils.get_index_at_wavelength(fiber_params('center_wavelength_nm'));
				n_rho = n_clad + fiber_params('nr_offset_from_cladding');
				n_rho_init = n_clad + init_fiber_params('nr_offset_from_cladding');
				n_change = n_rho - n_rho_init;

				subplot(231);
				plot(rho_arr*1e6, n_change*1e3, 'linewidth', 2); axis square; axis tight;
	            xlabel('r (\mum)'); ylabel('(n_{opt}-n_{init}) \times 10^{-3}');

				subplot(232);
				plot(rho_arr*1e6, n_rho, rho_arr*1e6, n_rho_init, 'linewidth', 2); axis square; axis tight;
	            xlabel('r (\mum)'); ylabel('n(r)');
	            legend('Optimized', 'Initial');

	            subplot(233);
				plot((fiber_params('neff') - n_clad)*1e3,'-+', 'linewidth', 2);
				xlabel('Mode index'); ylabel('(n_{eff} - n_{clad}) \times 10^{-3}');
				axis square; axis tight;

				utils.plot_cell_array('JLT figure 1', 2, 3, [4 5], 1:nn+1, MD_coeffs_hist, 'Iteration', 'Group delay (ps/m)');
				dsfig('JLT figure 1'); subplot(2,3,6);
				plot(1:nn+1, rms_MD_hist, 'linewidth', 2);
				xlabel('Iteration'); ylabel('rms group delay (ps/m)'); axis tight;
				
				%utils.plot_results(fiber_params, init_fiber_params);
				drawnow;
			end

			lambda_arr = linspace(1500,1600,20)*1e-9;
			rms_gd_wavelength_variation = test_robustness.vary_wavelength_rms_gd(fiber_params, lambda_arr);
			init_rms_gd_wavelength_variation = test_robustness.vary_wavelength_rms_gd(init_fiber_params, lambda_arr);
			dsfig('JLT figure 2');
			plot(lambda_arr*1e9, rms_gd_wavelength_variation, lambda_arr*1e9, init_rms_gd_wavelength_variation, 'linewidth', 2);
			xlabel('Wavelength (nm)'); ylabel('rms group delay (ps/m)')
			legend('Optimized', 'Initial');
			axis square; axis tight;
		end

		function obj = MMF_GI_freeform_reduce_coupling(obj)
			fiber_params = obj.bank.get_GI_MMF_1();
			fiber_params = utils.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 80;
			opt_params('max_dn') = 1e-4;
			opt_params('direction') =  "MIN"; %"MAX"; %"MIN";

			neff_hist = {fiber_params('neff')};

			for nn = 1:opt_params('opt_steps')
				disp(sprintf('Running iteration %d of %d...', nn, opt_params('opt_steps')));
				% compute and apply the index update
				index_update = obj.optimizer.opt_mode_coupling_freeform(fiber_params, init_fiber_params, opt_params);
				fiber_params('nxy_offset_from_cladding') = fiber_params('nxy_offset_from_cladding') + index_update;

				% update the realized fiber properties
				fiber_params = utils.solve_fiber_properties(fiber_params);

				% update the history
				neff_hist{end+1} = fiber_params('neff');
				utils.plot_cell_array('neff evolution', 1:nn+1, neff_hist, 'Iteration', 'Effective Index');
				utils.plot_results(fiber_params, init_fiber_params);

				drawnow;
			end

			lambda_arr = linspace(1200,1700,20)*1e-9;
			neff_wavelength_variation = test_robustness.vary_wavelength_neff(fiber_params, lambda_arr);
			utils.plot_cell_array('wavelength variation', lambda_arr*1e9, neff_wavelength_variation, 'Wavelength (nm)', 'Effective Index');
			axis tight;
		end

		function obj = MMF_GI_radial_reduce_coupling(obj)
			fiber_params = obj.bank.get_GI_MMF_2();
			fiber_params = utils.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 100;
			opt_params('max_dn') = 1e-4;
			opt_params('direction') =  "MIN"; %"MAX"; %"MIN";

			neff_hist = {fiber_params('neff')};

			for nn = 1:opt_params('opt_steps')
				disp(sprintf('Running iteration %d of %d...', nn, opt_params('opt_steps')));
				% compute and apply the index update
				index_update = obj.optimizer.opt_mode_coupling_radial(fiber_params, init_fiber_params, opt_params);
				fiber_params('nr_offset_from_cladding') = fiber_params('nr_offset_from_cladding') + index_update;

				% update the realized fiber properties
				fiber_params = utils.solve_fiber_properties(fiber_params);

				% update the history
				neff_hist{end+1} = fiber_params('neff');
				utils.plot_cell_array('neff evolution', 1,1,1, 1:nn+1, neff_hist, 'Iteration', 'Effective Index');
				utils.plot_results(fiber_params, init_fiber_params);

				drawnow;
			end
		end

		function obj = MMF_GI_radial_increase_coupling(obj)
			fiber_params = obj.bank.get_GI_MMF_1();
			fiber_params = utils.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 80;
			opt_params('max_dn') = 1e-4;
			opt_params('direction') =  "MAX"; 

			neff_hist = {fiber_params('neff')};

			for nn = 1:opt_params('opt_steps')
				disp(sprintf('Running iteration %d of %d...', nn, opt_params('opt_steps')));
				% compute and apply the index update
				index_update = obj.optimizer.opt_mode_coupling_freeform(fiber_params, init_fiber_params, opt_params);
				fiber_params('nxy_offset_from_cladding') = fiber_params('nxy_offset_from_cladding') + index_update;

				% update the realized fiber properties
				fiber_params = utils.solve_fiber_properties(fiber_params);

				% update the history
				neff_hist{end+1} = fiber_params('neff');
				utils.plot_cell_array('neff evolution', 1:nn+1, neff_hist, 'Iteration', 'Effective Index');
				utils.plot_results(fiber_params, init_fiber_params);

				drawnow;
			end

			lambda_arr = linspace(1200,1700,20)*1e-9;
			neff_wavelength_variation = test_robustness.vary_wavelength_neff(fiber_params, lambda_arr);
			utils.plot_cell_array('wavelength variation', lambda_arr*1e9, neff_wavelength_variation, 'Wavelength (nm)', 'Effective Index');
			axis tight;
		end

		%%%%%%%%%%%%%%%%%%%%%%%%%
		%   RESULTS NOT GREAT   %
		%%%%%%%%%%%%%%%%%%%%%%%%%

		function obj = MMF_SI_radial_reduce_coupling(obj)
			fiber_params = obj.bank.get_StepIndex_MMF();
			fiber_params = utils.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 1; %50;
			opt_params('max_dn') = 0*1e-4;
			opt_params('direction') =  "MIN"; %"MAX"; %"MIN";

			neff_hist = {fiber_params('neff')};

			for nn = 1:opt_params('opt_steps')
				disp(sprintf('Running iteration %d of %d...', nn, opt_params('opt_steps')));
				% compute and apply the index update
				index_update = obj.optimizer.opt_mode_coupling_radial(fiber_params, init_fiber_params, opt_params);
				fiber_params('nr_offset_from_cladding') = fiber_params('nr_offset_from_cladding') + index_update;

				% update the realized fiber properties
				fiber_params = utils.solve_fiber_properties(fiber_params);

				% update the history
				neff_hist{end+1} = fiber_params('neff');
				utils.plot_cell_array('neff evolution', 1:nn+1, neff_hist, 'Iteration', 'Effective Index');
				utils.plot_results(fiber_params, init_fiber_params);

				drawnow;
			end
		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%         NOT WORKING       %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		function obj = RCF_radial_reduce_coupling(obj) % NOT WORKING
			fiber_params = obj.bank.get_RCF_1();
			fiber_params = utils.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 100;
			opt_params('max_dn') = 0.5e-4;
			opt_params('direction') =  "RCF_v1";%"MIN";

			neff_hist = {fiber_params('neff')};

			for nn = 1:opt_params('opt_steps')
				disp(sprintf('Running iteration %d of %d...', nn, opt_params('opt_steps')));
				% compute and apply the index update
				index_update = obj.optimizer.opt_mode_coupling_radial(fiber_params, init_fiber_params, opt_params);
				fiber_params('nr_offset_from_cladding') = fiber_params('nr_offset_from_cladding') + index_update;

				% update the realized fiber properties
				fiber_params = utils.solve_fiber_properties(fiber_params);

				% update the history
				neff_hist{end+1} = fiber_params('neff');
				utils.plot_cell_array('neff evolution', 1, 1, 1, 1:nn+1, neff_hist, 'Iteration', 'n_{eff}');
				utils.plot_results(fiber_params, init_fiber_params);

				drawnow;
			end
		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%     UNDER DEVELOPMENT       %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		function obj = SMF_SI_increase_CD(obj)
			fiber_params = obj.bank.get_Corning_SMF28(); %obj.bank.get_GI_SMF_1();
			fiber_params = utils.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 50;
			opt_params('max_dn') = 0.5e-4;
			opt_params('direction') =  "MAX"; 

			rms_CD_hist = [rms(fiber_params('CD_coeffs_psnmkm'))];
			CD_coeffs_hist = {fiber_params('CD_coeffs_psnmkm')};

			for nn = 1:opt_params('opt_steps')
				disp(sprintf('Running iteration %d of %d...', nn, opt_params('opt_steps')));
				% compute and apply the index update
				index_update = obj.optimizer.opt_chromatic_dispersion_radial(fiber_params, init_fiber_params, opt_params);
				fiber_params('nr_offset_from_cladding') = fiber_params('nr_offset_from_cladding') + index_update;

				% update the realized fiber properties
				fiber_params = utils.solve_fiber_properties(fiber_params);

				% update the history
				rms_CD_hist = [rms_CD_hist rms(fiber_params('CD_coeffs_psnmkm'))];
				CD_coeffs_hist{end+1} = fiber_params('CD_coeffs_psnmkm');

				utils.plot_cell_array('CD evolution', 1:nn+1, CD_coeffs_hist, 'Iteration', 'Chromatic dispersion (ps/nm*km)');

				dsfig('rms CD evolution');
				plot(1:nn+1, rms_CD_hist);
				xlabel('Iteration'); ylabel('rms chromatic dispersion (ps/nm*km)'); axis tight;
				
				utils.plot_results(fiber_params, init_fiber_params);
				drawnow;
			end
		end

		function obj = run_playground(obj) % workspace for debugging and/or developing new fiber designs
			fiber_params = obj.bank.get_4C_MCF_debug();
			fiber_params = utils.solve_fiber_properties(fiber_params);

			utils.visualize_fiber(fiber_params);
		end

	end
end