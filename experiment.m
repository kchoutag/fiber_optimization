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

		function results = MingJunLi_vs_Karthik_1(obj)
			% run Ming Jung Li's fiber
			lambda_nm = 850;

			B = [];
            fib_orig = obj.bank.get_MJL_Alpha2p095(lambda_nm, B);
            rho_arr = linspace(0, fib_orig('nr')*fib_orig('dr'), fib_orig('nr'));


            B = [-0.0007 0.0007]; 
            %B = [0.0007]; 
            fib_pert = obj.bank.get_MJL_Alpha2p095(lambda_nm, B);

            dsfig('Index profiles');
 			subplot(121);
            plot(rho_arr*1e6, fib_orig('nr_offset_from_cladding'), ...
            	 rho_arr*1e6, fib_pert('nr_offset_from_cladding'), 'linewidth', 2);
            xlabel('Radial position r [\mum]'); ylabel('n(r)-n_{clad} (r)');
            legend('Pure alpha: n_{orig}(r)', 'Alpha + modified by B: n_{pert}(r)');
            grid on;

            subplot(122);
            plot(rho_arr*1e6, fib_pert('nr_offset_from_cladding') - fib_orig('nr_offset_from_cladding'), ...
            	'linewidth', 2);
             xlabel('Radial position [\mum]'); ylabel('n_{pert}(r)-n_{orig} (r)');
             grid on;
            
            dsfig('MJL Patent Fig. 6C');
            plot(rho_arr*1e6, fib_pert('nr_non_alpha_part'), 'linewidth', 2);
            xlabel('Radial position [\mum]'); ylabel('H(x) = G(x)-x^{a0/2}');
            grid on;

			D_valid = 120; % MJL's patent says only first 15 mode groups are considered. 15*16/2 = 120 modes total.
			
            results = [];
            
            % compare original and perturbed fiber
            fib_orig = utils.solve_fiber_properties(fib_orig);
			gds = fib_orig('MD_coeffs_psm');
			gd_orig = rms(gds(1:D_valid))

            
			fib_pert = utils.solve_fiber_properties(fib_pert);
			gds = fib_pert('MD_coeffs_psm');
			gd_pert = rms(gds(1:D_valid))

			utils.MJL_visualize_fiber(fib_pert);
		end

		function obj = MingJunLi_vs_Karthik_2(obj)
			B = [-0.0007 0.0007]; 
			lambda_nm = 850;

			fiber_params = obj.bank.get_MJL_Alpha2p095(lambda_nm, B);
			fiber_params = utils.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 100;
			opt_params('max_dn') = 1e-5;
			opt_params('direction') =  "MIN"; 
			opt_params('apply_smoothing') = true;
			opt_params('gaussian_filter_std_m') = 1e-6;


			rms_MD_hist = [rms(fiber_params('MD_coeffs_psm'))];
			MD_coeffs_hist = {fiber_params('MD_coeffs_psm')};

			rho_arr = linspace(0, fiber_params('nr')*fiber_params('dr'), fiber_params('nr'));
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
				dsfig('MJL fiber with my design method');
	        	n_clad = utils.get_index_at_wavelength(fiber_params('center_wavelength_nm'));
				n_rho = n_clad + fiber_params('nr_offset_from_cladding');
				n_rho_init = n_clad + init_fiber_params('nr_offset_from_cladding');
				n_change = n_rho - n_rho_init;

				subplot(231);
				plot(rho_arr*1e6, n_change, 'k', 'linewidth', 2); axis square; axis tight;
	            xlabel('r (\mum)'); ylabel('n_{opt}-n_{init}');

				subplot(232);
				plot(rho_arr*1e6, n_rho, rho_arr*1e6, n_rho_init, 'linewidth', 2); axis square; axis tight;
	            xlabel('r (\mum)'); ylabel('n(r)');
	            legend('Optimized', 'Initial');

	            subplot(233);
				plot((fiber_params('neff') - n_clad),'k-+', 'linewidth', 2);
				xlabel('Mode index'); ylabel('n_{eff} - n_{clad}');
				ylim([0 max(fiber_params('neff') - n_clad)]);
				xlim([1 fiber_params('D')]);
				axis square; 

				utils.plot_cell_array('MJL fiber with my design method', 2, 3, [4 5], 1:nn+1, MD_coeffs_hist, 'Iteration', 'Group delay (ps/m)', 'line');
				dsfig('MJL fiber with my design method'); subplot(2,3,6);
				plot(0:nn, rms_MD_hist, 'k', 'linewidth', 2);
				xlabel('Iteration'); ylabel('rms group delay (ps/m)'); axis tight;
				ylim([0 max(rms_MD_hist)]);
				
				%utils.plot_results(fiber_params, init_fiber_params);
				drawnow;
			end

			%save the results
			if(1)
				save('saved_fibers/MJL_Karthik_optimized_fiber.mat');
			end
		end

		function obj = test_robustness_GI_MMF_low_MD(obj)
			load('saved_fibers/GI_MMF_low_MD.mat');


            % Used this plot to figure out the max spatial frequency is ~0.25e6 1/m
            % figure;
            % fx = 1/(rho_arr(2) - rho_arr(1));
            % plot(linspace(-fx/2, fx/2, length(n_change)), abs(fftshift(fft(n_change))));

        	rng(101);

            dsfig('Robustness: low frequency perturbations');
            opt_profile = fiber_params('nr_offset_from_cladding');
            opt_change = n_change;
            
            subplot(131);
            plot(1:10);
            axis square;
            title('Placeholder for f rescaling plot');

            % EXAMPLE OF PERTURBATION
            subplot(132); 
            plot(rho_arr*1e6, opt_change, 'linewidth', 2, 'DisplayName', sprintf('Optimal, rms GD = %.2f ps/km', 1e3*rms(fiber_params('MD_coeffs_psm')))); 
			hold on;
			n_terms = 5;
			f_pert = linspace(0.1, 0.25e6, n_terms); % used the commented out plot above to figure out max spatial frequency is ~0.25e6 1/m
			n_pert = zeros(size(opt_change));
			pert_amplitude = 2e-5;
			for ii = 1:n_terms
				n_pert = n_pert + ...
				(randn()*cos(2*pi*f_pert(ii)*rho_arr) + randn()*sin(2*pi*f_pert(ii)*rho_arr));
            end
            % normalize perturbation to optimal profile to zero-mean with amplitude pert_amplitude
            n_pert = n_pert - mean(n_pert);
			n_pert = pert_amplitude * (n_pert / rms(n_pert));
			
			fiber_params('nr_offset_from_cladding') = opt_profile + n_pert;
			fiber_params = utils.solve_fiber_properties(fiber_params);

			plot(rho_arr*1e6, opt_change + n_pert, 'linewidth', 2, 'DisplayName', sprintf('Optimal + random deviations, rms GD = %.2f ps/km', 1e3*rms(fiber_params('MD_coeffs_psm')))); 

			axis square; axis tight;
            xlabel('r (\mum)'); ylabel('n_{opt}-n_{init}');
            legend;

            % SWEEP PERTURBATION STRENGTHS
            subplot(133); 
            MAX_PERT_STD = 5e-5;
            n_fibers_per_pert = 20;
            n_pert_strengths = 10;
            pert_amps = MAX_PERT_STD * (1/n_pert_strengths : 1/n_pert_strengths : 1);
            rms_gd_table = zeros(n_fibers_per_pert, n_pert_strengths);

            for ii = 1:n_pert_strengths
            	fprintf('Testing perturbation strength #%d of %d...\n', ii, n_pert_strengths);
            	for nn = 1:n_fibers_per_pert
            		fprintf(' ...fiber #%d of %d\n', nn, n_fibers_per_pert);
		            n_pert = zeros(size(n_change));
					pert_amplitude = pert_amps(ii);
					for kk = 1:n_terms
						n_pert = n_pert + ...
						(randn()*cos(2*pi*f_pert(kk)*rho_arr) + randn()*sin(2*pi*f_pert(kk)*rho_arr));
					end
					% normalize perturbation to optimal profile to zero-mean with amplitude pert_amplitude
		            n_pert = n_pert - mean(n_pert);
					n_pert = pert_amplitude * (n_pert / rms(n_pert));
					
					fiber_params('nr_offset_from_cladding') = opt_profile + n_pert;
					fiber_params = utils.solve_fiber_properties(fiber_params);

					rms_gd_table(nn, ii) = rms(fiber_params('MD_coeffs_psm'));
				end
			end
            %boxplot(1e3*rms_gd_table, pert_amps);
            plot(pert_amps, 1e3*mean(rms_gd_table,1), '-ko', 'linewidth', 2);
			xlabel('rms amplitude of random deviation from optimal profile');
			ylabel('Average rms GD of MMF (ps/km)');
			axis square; axis tight;
			ylims = ylim;
			ylim([0 ylims(2)]);

			% save the results
			if(0)
				save('saved_fibers/GI_MMF_low_MD_robustness.mat');
			end
		end

		function obj = MMF_GI_radial_reduce_modal_dispersion(obj)
			fiber_params = obj.bank.get_MMF_withTrench_1();
			fiber_params = utils.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 200;
			opt_params('max_dn') = 1e-5;
			opt_params('direction') =  "MIN"; 
			opt_params('apply_smoothing') = true;
			opt_params('gaussian_filter_std_m') = 1e-6;


			rms_MD_hist = [rms(fiber_params('MD_coeffs_psm'))];
			MD_coeffs_hist = {fiber_params('MD_coeffs_psm')};

			rho_arr = linspace(0, fiber_params('nr')*fiber_params('dr'), fiber_params('nr'));
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
	        	n_clad = utils.get_index_at_wavelength(fiber_params('center_wavelength_nm'));
				n_rho = n_clad + fiber_params('nr_offset_from_cladding');
				n_rho_init = n_clad + init_fiber_params('nr_offset_from_cladding');
				n_change = n_rho - n_rho_init;

				subplot(231);
				plot(rho_arr*1e6, n_change, 'k', 'linewidth', 2); axis square; axis tight;
	            xlabel('r (\mum)'); ylabel('n_{opt}-n_{init}');

				subplot(232);
				plot(rho_arr*1e6, n_rho, rho_arr*1e6, n_rho_init, 'linewidth', 2); axis square; axis tight;
	            xlabel('r (\mum)'); ylabel('n(r)');
	            legend('Optimized', 'Initial');

	            subplot(233);
				plot((fiber_params('neff') - n_clad),'k-+', 'linewidth', 2);
				xlabel('Mode index'); ylabel('n_{eff} - n_{clad}');
				ylim([0 max(fiber_params('neff') - n_clad)]);
				xlim([1 fiber_params('D')]);
				axis square; 

				utils.plot_cell_array('JLT figure 1', 2, 3, [4 5], 1:nn+1, MD_coeffs_hist, 'Iteration', 'Group delay (ps/m)', 'line');
				dsfig('JLT figure 1'); subplot(2,3,6);
				plot(0:nn, rms_MD_hist, 'k', 'linewidth', 2);
				xlabel('Iteration'); ylabel('rms group delay (ps/m)'); axis tight;
				ylim([0 max(rms_MD_hist)]);
				
				%utils.plot_results(fiber_params, init_fiber_params);
				drawnow;
			end
			interp_factor = 1.5;
			fiber_bigger_grid = test_robustness.change_gridding(fiber_params, interp_factor);
			fiber_bigger_grid = utils.solve_fiber_properties(fiber_bigger_grid);
			dsfig('Changing Grid');
			subplot(121);
			plot(fiber_bigger_grid('nr_offset_from_cladding'));
			title('Shape of index profile in new grid');
			subplot(122);
			bar(fiber_bigger_grid('MD_coeffs_psm'));
			title('Group delays (ps/m)')
			drawnow;

			lambda_arr = linspace(1400,1900,30)*1e-9;
			rms_gd_wavelength_variation = test_robustness.vary_wavelength_rms_gd(fiber_params, lambda_arr);
			init_rms_gd_wavelength_variation = test_robustness.vary_wavelength_rms_gd(init_fiber_params, lambda_arr);
			dsfig('JLT figure 2');
			subplot(121);
			plot(lambda_arr*1e9, rms_gd_wavelength_variation, lambda_arr*1e9, init_rms_gd_wavelength_variation, 'linewidth', 2);
			xlabel('Wavelength (nm)'); ylabel('rms group delay (ps/m)')
			legend('Optimized', 'Initial');
			axis square; axis tight;
			D_wavelength_variation = test_robustness.vary_wavelength_number_modes(fiber_params, lambda_arr);
			init_D_wavelength_variation = test_robustness.vary_wavelength_number_modes(init_fiber_params, lambda_arr);
			dsfig('JLT figure 2');
			subplot(122);
			plot(lambda_arr*1e9, D_wavelength_variation, lambda_arr*1e9, init_D_wavelength_variation, 'linewidth', 2);
			xlabel('Wavelength (nm)'); ylabel('Number of modes D')
			legend('Optimized', 'Initial');
			axis square; axis tight; %ylim([0, max(max(D_wavelength_variation), max(init_D_wavelength_variation))+1]);

			%save the results
			if(1)
				save('saved_fibers/GI_MMF_low_MD.mat');
			end
		end

		function obj = MMF_GI_radial_reduce_modal_dispersion_part2(obj)
			% continue analysis of optimization without designing the fiber from scratch again
			load('saved_fibers/GI_MMF_low_MD.mat');

			figure();

			x = init_fiber_params('nr_offset_from_cladding');
			y = fiber_params('nr_offset_from_cladding');
 
			diff_nr = y - x;
			plot(diff_nr);

			sweep_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy
			rms_gd_variation = [];
			r_arr =  0:0.01:1.5;
			for rr = 1:length(r_arr)
				r = r_arr(rr);
				disp(sprintf('testing sensitivity of optimization %d of %d...', rr, length(r_arr)));
				sweep_fiber_params('nr_offset_from_cladding') = init_fiber_params('nr_offset_from_cladding') + r*diff_nr;
				sweep_fiber_params = utils.solve_fiber_properties(sweep_fiber_params);
				rms_gd_variation = [rms_gd_variation rms(sweep_fiber_params('MD_coeffs_psm'))];
			end

			plot(r_arr, rms_gd_variation, 'k', 'linewidth', 2);
			xlabel('Deviation parameter f'); ylabel('rms group delay (ps/m)');
			axis tight; axis square;
			ylim([0 max(rms_gd_variation)]);
		end

		function obj = RCF_freeform_increase_degeneracies(obj) % see if RCF -> MCF when we try to increase the modal degeneracies
			fiber_params = obj.bank.get_RCF_3(); %obj.bank.get_RCF_2();
			fiber_params = utils.solve_fiber_properties(fiber_params);

			init_fiber_params = containers.Map(fiber_params.keys, fiber_params.values); % make copy of the init fiber params container

			opt_params = containers.Map;
			opt_params('opt_steps') = 380;
			opt_params('max_dn') = 2e-4;
			opt_params('direction') =  "RCF_to_MCF";

			neff_hist = {fiber_params('neff') - fiber_params('n_clad')};

			for nn = 1:opt_params('opt_steps')
				disp(sprintf('Running iteration %d of %d...', nn, opt_params('opt_steps')));
				% compute and apply the index update
				index_update = obj.optimizer.opt_mode_coupling_freeform(fiber_params, init_fiber_params, opt_params);

				fiber_params('nxy_offset_from_cladding') = fiber_params('nxy_offset_from_cladding') + index_update;

				% update the realized fiber properties
				fiber_params = utils.solve_fiber_properties(fiber_params);

				% update the history
				neff_hist{end+1} = fiber_params('neff')- fiber_params('n_clad');
				utils.plot_cell_array('neff evolution for JLT', 1, 1, 1, 1:nn+1, neff_hist , 'Iteration', 'n_{eff}-n_{clad}', 'line');
				pbaspect([6 2 1]);
				utils.plot_results(fiber_params, init_fiber_params);

				drawnow;
			end

			% plot the initial and final modal fields
			dsfig('modal fields for JLT');
			init_fields = init_fiber_params('fields'); opt_fields = fiber_params('fields');
			D_init = init_fiber_params('D'); D_opt = fiber_params('D');
			[ny, nx] = size(opt_fields(:,:,1));
			allfields = zeros(2*nx, max(D_init,D_opt)*ny);
			for ii = 1:D_init
				nr = 1;
				nc = ii;
				allfields((nr-1)*nx+1:nr*nx, (nc-1)*nx+1 : nc*nx) = abs(init_fields(:,:,ii)).^2;
			end
			for ii = 1:D_opt
				nr = 2;
				nc = ii;
				allfields((nr-1)*nx+1:nr*nx, (nc-1)*nx+1 : nc*nx) = abs(opt_fields(:,:,ii)).^2;
			end
			imagesc(allfields);
			pbaspect([max(D_init,D_opt) 2 1]); colorbar; axis off;
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
				utils.plot_cell_array('neff evolution', 1:nn+1, neff_hist, 'Iteration', 'Effective Index', 'line');
				utils.plot_results(fiber_params, init_fiber_params);

				drawnow;
			end

			lambda_arr = linspace(1200,1700,20)*1e-9;
			neff_wavelength_variation = test_robustness.vary_wavelength_neff(fiber_params, lambda_arr);
			utils.plot_cell_array('wavelength variation', lambda_arr*1e9, neff_wavelength_variation, 'Wavelength (nm)', 'Effective Index', 'line');
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
				utils.plot_cell_array('neff evolution', 1,1,1, 1:nn+1, neff_hist, 'Iteration', 'Effective Index', 'line');
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
				utils.plot_cell_array('neff evolution', 1:nn+1, neff_hist, 'Iteration', 'Effective Index', 'line');
				utils.plot_results(fiber_params, init_fiber_params);

				drawnow;
			end

			lambda_arr = linspace(1200,1700,20)*1e-9;
			neff_wavelength_variation = test_robustness.vary_wavelength_neff(fiber_params, lambda_arr);
			utils.plot_cell_array('wavelength variation', lambda_arr*1e9, neff_wavelength_variation, 'Wavelength (nm)', 'Effective Index', 'line');
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
				utils.plot_cell_array('neff evolution', 1:nn+1, neff_hist, 'Iteration', 'Effective Index', 'line');
				utils.plot_results(fiber_params, init_fiber_params);

				drawnow;
			end
		end

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%         NOT WORKING       %
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		

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

				utils.plot_cell_array('CD evolution', 1,1,1, 1:nn+1, CD_coeffs_hist, 'Iteration', 'Chromatic dispersion (ps/nm*km)', 'line');

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