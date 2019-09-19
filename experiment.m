classdef experiment
	properties
		bank = fiber_bank();
		
		optimizer = optimizers();
	end
	methods
		function obj = experiment()
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
			opt_params('opt_steps') = 50;
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
				utils.plot_cell_array('neff evolution', 1:nn+1, neff_hist, 'Iteration', 'Effective Index');
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

		function obj = MMF_GI_radial_reduce_modal_dispersion(obj)
			fiber_params = obj.bank.get_MMF_withTrench_1();
			fiber_params = utils.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 80;
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

				utils.plot_cell_array('MD evolution', 1:nn+1, MD_coeffs_hist, 'Iteration', 'Group delay (ps/m)');

				dsfig('rms MD evolution');
				plot(1:nn+1, rms_MD_hist, 'linewidth', 2);
				xlabel('Iteration'); ylabel('rms group delay (ps/m)'); axis tight;
				
				utils.plot_results(fiber_params, init_fiber_params);
				drawnow;
			end

			lambda_arr = linspace(1500,1600,20)*1e-9;
			rms_gd_wavelength_variation = test_robustness.vary_wavelength_rms_gd(fiber_params, lambda_arr);
			utils.plot_cell_array('Wavelength Variation', lambda_arr*1e9, rms_gd_wavelength_variation, 'Wavelength (nm)', 'rms group delay (ps/m)');
			axis tight;
		end

		function obj = SMF_SI_increase_CD()
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
			fiber_params = obj.bank.get_GI_MMF_2();
			fiber_params = utils.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 80;
			opt_params('max_dn') = 0.5e-4;
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

				utils.plot_cell_array('MD evolution', 1:nn+1, MD_coeffs_hist, 'Iteration', 'Group delay (ps/m)');

				dsfig('rms MD evolution');
				plot(1:nn+1, rms_MD_hist);
				xlabel('Iteration'); ylabel('rms group delay (ps/m)'); axis tight;
				
				utils.plot_results(fiber_params, init_fiber_params);
				drawnow;
			end
		end

	end
end