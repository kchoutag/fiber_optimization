classdef utils
	methods(Static)
		function plot_cell_array(fig_name, fig_rows, fig_cols, subplot_locs, x_vals, cell_in, x_label, y_label, style)
			dsfig(fig_name);
			subplot(fig_rows, fig_cols, subplot_locs);

			dim1 = length(cell_in);
			dim2 = 0;
			for ii = 1:dim1
				dim2 = max(dim2, length(cell_in{ii}));
			end
			cell_array = nan(dim1,dim2);
			for ii = 1:dim1
				cell_array(ii,1:length(cell_in{ii})) = sort(cell_in{ii},'descend');
            end
            
            if(length(x_vals) > 1)
            	if(strcmp(style,'line'))
                	plot(x_vals, cell_array,'-k','linewidth',2);
                elseif(strcmp(style, 'circle'))
                	plot(x_vals, cell_array,'ko','markersize',3,'markerfacecolor', 'k');
                elseif(strcmp(style, 'plus'))
                	plot(x_vals, cell_array,'+','markersize',2);
                end
                xlabel(x_label);
                ylabel(y_label);
                axis tight;
            end
		end

		function visualize_fiber(fiber_params)
			dsfig('Fiber index profile');

			if(~fiber_params('axially_symm'))
				x_arr = linspace(-fiber_params('dx')*fiber_params('nx')*1e6/2, fiber_params('dx')*fiber_params('nx')*1e6/2, fiber_params('nx'));
				y_arr = linspace(-fiber_params('dy')*fiber_params('ny')*1e6/2, fiber_params('dy')*fiber_params('ny')*1e6/2, fiber_params('ny'));
				
				subplot(121);
				n_clad = utils.get_index_at_wavelength(fiber_params('center_wavelength_nm'));
				n_xy = n_clad + fiber_params('nxy_offset_from_cladding');
				imagesc(x_arr, y_arr, n_xy); axis square; colorbar;
	            xlabel('x (\mum)'); ylabel(' y (\mum)');

	            subplot(122);
	            [X, Y] = meshgrid(x_arr, y_arr);
				surf(X, Y, n_xy); colorbar;
				axis square; shading interp;
				xlabel('x (\mum)'); ylabel('y (\mum)');
	        else 
	        	rho_arr = linspace(0, fiber_params('nr')*fiber_params('dr'), fiber_params('nr'));
	        	n_clad = utils.get_index_at_wavelength(fiber_params('center_wavelength_nm'));
				n_rho = n_clad + fiber_params('nr_offset_from_cladding');

				subplot(121);
				plot(rho_arr*1e6, (n_rho-n_clad)*1e3, 'linewidth', 2); axis square;
	            xlabel('r (\mum)'); ylabel('n(r) - n_{clad} \times 10^{-3}');

	            subplot(122);
	            [n_xy, x_arr, y_arr] = profiles.synthesize_nxy_from_nr(n_rho, rho_arr, n_clad);
	            [X, Y] = meshgrid(x_arr, y_arr);
	            surf(X*1e6, Y*1e6, n_xy); colorbar;
				axis square; shading interp;
				xlabel('x (\mum)'); ylabel('y (\mum)');
	        end
	        	
	       	dsfig('Fiber properties');
	       	
	       	subplot(221);
			plot(fiber_params('Aeff'), '-+', 'linewidth', 2);
			xlabel('Mode index'); ylabel('Effective area (\mum^2)');
			axis square; axis tight;

			subplot(222);
			plot(fiber_params('MD_coeffs_psm'), '-+', 'linewidth', 2);
			xlabel('Mode index'); ylabel('Group delays (ps/m)');
			axis square; axis tight;

			subplot(223);
			plot(fiber_params('CD_coeffs_psnmkm'), '-+', 'linewidth', 2);
			xlabel('Mode index'); ylabel('Chromatic dispersion (ps/nm*km)');
			axis square; axis tight;

			subplot(224);
			plot((fiber_params('neff') - n_clad)*1e3, '-+', 'linewidth', 2);
			xlabel('Mode index'); ylabel('n_{eff} - n_{clad} (\times 10^{-3})');
			axis square; axis tight;


			dsfig('Modes of Fiber');
			nx_modes = floor(sqrt(fiber_params('D')));
			ny_modes = ceil(fiber_params('D')/nx_modes);
			opt_fields = fiber_params('fields'); 
			[ny, nx] = size(opt_fields(:,:,1));

			allfields = zeros(nx_modes*nx, ny_modes*ny);
			for ii = 1:fiber_params('D')
				nr = floor((ii-1)/ny_modes) + 1;
				nc = mod(ii-1,ny_modes)+1;
				allfields((nr-1)*nx+1:nr*nx, (nc-1)*nx+1 : nc*nx) = abs(opt_fields(:,:,ii)).^2;
			end
			imagesc(allfields);
			pbaspect([ny_modes nx_modes 1]); colorbar; axis off;
		end

		function plot_results(opt_fiber_params, init_fiber_params)
			dsfig('Index profile');

			if(~opt_fiber_params('axially_symm'))
				x_arr = linspace(-opt_fiber_params('dx')*opt_fiber_params('nx')*1e6/2, opt_fiber_params('dx')*opt_fiber_params('nx')*1e6/2, opt_fiber_params('nx'));
				y_arr = linspace(-opt_fiber_params('dy')*opt_fiber_params('ny')*1e6/2, opt_fiber_params('dy')*opt_fiber_params('ny')*1e6/2, opt_fiber_params('ny'));
				n_clad = utils.get_index_at_wavelength(opt_fiber_params('center_wavelength_nm'));

				subplot(131);
				n_init = n_clad + init_fiber_params('nxy_offset_from_cladding');
				imagesc(x_arr, y_arr, n_init); axis square; colorbar;
	            xlabel('x (\mum)'); ylabel(' y (\mum)');
	            title('n_{init}');

				subplot(132);
				n_change = opt_fiber_params('nxy_offset_from_cladding') - init_fiber_params('nxy_offset_from_cladding');
				imagesc(x_arr, y_arr, n_change); axis square; colorbar;
	            xlabel('x (\mum)'); ylabel(' y (\mum)');
	            title('n_{opt} - n_{init}');

				subplot(133);
				n_xy = n_clad + opt_fiber_params('nxy_offset_from_cladding');
				imagesc(x_arr, y_arr, n_xy); axis square; colorbar;
	            xlabel('x (\mum)'); ylabel(' y (\mum)');
	            title('n_{opt}');

	            %subplot(133);
	            %[X, Y] = meshgrid(x_arr, y_arr);
				%surf(X, Y, n_xy); colorbar;
				%axis square; shading interp;
				%xlabel('x (\mum)'); ylabel('y (\mum)');
	        else 
	        	rho_arr = linspace(0, opt_fiber_params('nr')*opt_fiber_params('dr'), opt_fiber_params('nr'));
	        	n_clad = utils.get_index_at_wavelength(opt_fiber_params('center_wavelength_nm'));
				n_rho = n_clad + opt_fiber_params('nr_offset_from_cladding');
				n_rho_init = n_clad + init_fiber_params('nr_offset_from_cladding');
				n_change = n_rho - n_rho_init;

				subplot(121);
				plot(rho_arr*1e6, n_change*1e3, 'linewidth', 2); axis square; axis tight;
	            xlabel('r (\mum)'); ylabel('n_{opt}-n_{init} \times 10^{-3}');

				subplot(122);
				plot(rho_arr*1e6, n_rho, rho_arr*1e6, n_rho_init, 'linewidth', 2); axis square; axis tight;
	            xlabel('r (\mum)'); ylabel('n(r)');
	            legend('Optimized', 'Initial');

	            %subplot(133);
	            %[n_xy, x_arr, y_arr] = profiles.synthesize_nxy_from_nr(n_rho, rho_arr, n_clad);
	            %[X, Y] = meshgrid(x_arr, y_arr);
	            %surf(X*1e6, Y*1e6, n_xy); colorbar;
				%axis square; shading interp;
				%xlabel('x (\mum)'); ylabel('y (\mum)');
	        end
	        	
	       	dsfig('Fiber Properties');
	       
	       	subplot(141);
			plot(opt_fiber_params('Aeff'), '-+', 'linewidth', 2);
			xlabel('Mode index'); ylabel('Effective area (\mum^2)');
			axis square; axis tight;

			subplot(142);
			plot(opt_fiber_params('MD_coeffs_psm'),'-+', 'linewidth', 2);
			xlabel('Mode index'); ylabel('Group delays (ps/m)');
			axis square; axis tight;

			subplot(143);
			plot(opt_fiber_params('CD_coeffs_psnmkm'),'-+', 'linewidth', 2);
			xlabel('Mode index'); ylabel('Chromatic dispersion (ps/nm*km)');
			axis square; axis tight;

			subplot(144);
			plot((opt_fiber_params('neff') - n_clad)*1e3,'-+', 'linewidth', 2);
			xlabel('Mode index'); ylabel('n_{eff} - n_{clad} (\times 10^{-3})');
			axis square; axis tight;

			dsfig('Modes');
			nx_modes = floor(sqrt(opt_fiber_params('D')));
			ny_modes = ceil(opt_fiber_params('D')/nx_modes);
			opt_fields = opt_fiber_params('fields'); 
			[ny, nx] = size(opt_fields(:,:,1));

			allfields = zeros(nx_modes*nx, ny_modes*ny);
			for ii = 1:opt_fiber_params('D')
				nr = floor((ii-1)/ny_modes) + 1;
				nc = mod(ii-1,ny_modes)+1;
				allfields((nr-1)*nx+1:nr*nx, (nc-1)*nx+1 : nc*nx) = abs(opt_fields(:,:,ii)).^2;
			end
			imagesc(allfields);
			pbaspect([ny_modes nx_modes 1]); colorbar; axis off;
		end

		function [n_fusedSi] = get_index_at_wavelength(wavelength_nm)
			% Sellmeier coefficients for fused Si - G. P. Agarwal
			B1 = 0.6961663;
			B2 = 0.4079426;
			B3 = 0.8974794;

			L1 = 0.0684043; %um
			L2 = 0.1162414; %um
			L3 = 9.896161;  %um

			c = utils.get_speed_light();

			omega1 = 2*pi*c/(L1*1e-6);
			omega2 = 2*pi*c/(L2*1e-6);
			omega3 = 2*pi*c/(L3*1e-6);
			omega  = 2*pi*c/(wavelength_nm*1e-9);
			n_fusedSi = sqrt(1 + B1*((omega1^2)/(omega1^2 - omega^2)) + ...
					        B2*((omega2^2)/(omega2^2 - omega^2)) + ...
					        B3*((omega3^2)/(omega3^2 - omega^2)));
		end

		function modal_intensity = integrate_mode_on_circle(modal_field, x_arr, y_arr, rho)
			% integrate |modal_field|^2 on a circle of radius rho
			tol = 0.5e-6;
			modal_intensity = 0;
			for xx = 1:length(x_arr)
                xval = x_arr(xx);
				for yy = 1:length(y_arr)
                    yval = y_arr(yy);
					rho_curr = sqrt(xval^2 + yval^2);
					if (abs(rho_curr - rho) <= tol)
						modal_intensity = modal_intensity + abs(modal_field(yy,xx))^2;
					end
				end
			end
		end

		function [fib_params] = solve_fiber_properties(fib_params)
			addpath('../modesolver-2011-04-22/');
			n_modes_upper_lim = 70;
			tol_neff = 0.15*1e-3; % to remove very-weakly guided modes
			c = utils.get_speed_light();
			lambda_nm = fib_params('center_wavelength_nm');
			d_lambda_nm = fib_params('d_wavelength_nm');

			if(~isKey(fib_params,'axially_symm'))
				error('Fiber parameters map does not contain axially_symm key!');
			end

			% center wavelength
			n_clad = utils.get_index_at_wavelength(lambda_nm);
			fib_params('n_clad') = n_clad;
			if(fib_params('axially_symm'))
				rho_arr = linspace(0, fib_params('nr')*fib_params('dr'), fib_params('nr'));
				[nxy_offset_from_cladding, x_arr, y_arr] = profiles.synthesize_nxy_from_nr(fib_params('nr_offset_from_cladding'), rho_arr, 0);
				dx = x_arr(2) - x_arr(1); dy = dx;
			else
				nxy_offset_from_cladding = fib_params('nxy_offset_from_cladding');
				dx = fib_params('dx'); dy = fib_params('dy');
			end
			nxy = nxy_offset_from_cladding + n_clad;
			n_core = max(max(nxy)); 
			eps_profile = (nxy).^2;
			[fields, neff_aux] = svmodes(lambda_nm*1e-9, n_core, n_modes_upper_lim, dx, dy, eps_profile, '0000', 'scalar');
			valid_idx = neff_aux > n_clad + tol_neff;
			fields = fields(:,:,valid_idx);
			neff = neff_aux(valid_idx);
			D_center = length(neff);

			% left wavelength
			left_wavelength_nm = lambda_nm - d_lambda_nm;
			n_clad = utils.get_index_at_wavelength(left_wavelength_nm);
			if(fib_params('axially_symm'))
				rho_arr = linspace(0, fib_params('nr')*fib_params('dr'), fib_params('nr'));
				[nxy_offset_from_cladding, x_arr, y_arr] = profiles.synthesize_nxy_from_nr(fib_params('nr_offset_from_cladding'), rho_arr, 0);
				dx = x_arr(2) - x_arr(1); dy = dx;
			else
				nxy_offset_from_cladding = fib_params('nxy_offset_from_cladding');
				dx = fib_params('dx'); dy = fib_params('dy');
			end
			nxy = nxy_offset_from_cladding + n_clad;
			n_core = max(max(nxy)); 
			eps_profile = (nxy).^2;
			[fields_left, neff_aux_left] = svmodes(left_wavelength_nm*1e-9, n_core, n_modes_upper_lim, dx, dy, eps_profile, '0000', 'scalar');
			valid_idx = neff_aux_left > n_clad + tol_neff;
			fields_left = fields_left(:,:,valid_idx);
			neff_left = neff_aux_left(valid_idx);
			D_left = length(neff_left);

			% right wavelength
			right_wavelength_nm = lambda_nm + d_lambda_nm;
			n_clad = utils.get_index_at_wavelength(right_wavelength_nm);
			if(fib_params('axially_symm'))
				rho_arr = linspace(0, fib_params('nr')*fib_params('dr'), fib_params('nr'));
				[nxy_offset_from_cladding, x_arr, y_arr] = profiles.synthesize_nxy_from_nr(fib_params('nr_offset_from_cladding'), rho_arr, 0);
				dx = x_arr(2) - x_arr(1); dy = dx;
			else
				nxy_offset_from_cladding = fib_params('nxy_offset_from_cladding');
				dx = fib_params('dx'); dy = fib_params('dy');
			end
			nxy = nxy_offset_from_cladding + n_clad;
			n_core = max(max(nxy)); 
			eps_profile = (nxy).^2;
			[fields_right, neff_aux_right] = svmodes(right_wavelength_nm*1e-9, n_core, n_modes_upper_lim, dx, dy, eps_profile, '0000', 'scalar');
			valid_idx = neff_aux_right > n_clad + tol_neff;
			fields_right = fields_right(:,:,valid_idx);
			neff_right = neff_aux_right(valid_idx);
			D_right = length(neff_right);

			if ((D_center ~= D_left) || (D_center ~= D_right) || (D_left ~= D_right))
				error('Number of modes not equal at left, center, right wavelengths!');
			end

			D = min(min(D_center, D_left), D_right);
			
			% calculate neff = beta0 of modes = phase velocities
			neff = neff(1:D);
			neff_left = neff_left(1:D);
			neff_right = neff_right(1:D);

			% set the fields
			fields = fields(:,:,1:D);
			fields_left = fields_left(:,:,1:D);
			fields_right = fields_right(:,:,1:D);

			% calculate effective areas (in [um^2] units)
			Aeff = zeros(1,D);
			for ii = 1:D
				% Aeff in terms of um^2
				Aeff(ii) = 1e6*1e6 * sum(sum(abs(fields(:,:,ii)).^2 * dx * dy))^2 / sum(sum(abs(fields(:,:,ii)).^4 * dx * dy));
			end
			fib_params('Aeff') = Aeff;
			
    		d_omega = (-2*pi*c/((lambda_nm*1e-9)^2)) * d_lambda_nm * 1e-9;
    		k0 = 2*pi/(lambda_nm*1e-9); k0_left = 2*pi/(left_wavelength_nm*1e-9); k0_right = 2*pi/(right_wavelength_nm*1e-9);
			beta_c = k0*neff; beta_left = k0_left*neff_left; beta_right = k0_right*neff_right;

			% calculate the modal dispersion for all modes using the central difference approximation (in [ps/m] units)
			modal_dispersion_coeffs_sm = (beta_left - beta_right)/(2*d_omega); %units s/m, note it's left-right since small omega corresponds to large lambda
			scale_sm_to_psm = 1e12;
			modal_dispersion_coeffs_psm = modal_dispersion_coeffs_sm * scale_sm_to_psm; % units ps/m
			modal_dispersion_coeffs_psm = modal_dispersion_coeffs_psm - mean(modal_dispersion_coeffs_psm); % make zero mean

    		% calculate chromatic dispersion for all modes (in [ps/nm*km] units)
			% using Sercan's code from /afs/ir/users/k/c/kchoutag/private/kahngroup/b/soarik/codes/GD compensation/dispersion.m
			scale_smm_to_psnmkm = 1e6; % conversion from s/(m*m) to ps/(nm*km)
			omega_arr = 2*pi*c ./ (1e-9*[right_wavelength_nm lambda_nm left_wavelength_nm]);
			chromatic_dispersion_psnmkm = zeros(1,D);
			for mm = 1:D
				beta_arr = [beta_right(mm) beta_c(mm) beta_left(mm)];
				[pdispersion,~,mudispersion]=polyfit(omega_arr,beta_arr,1);
				gd = pdispersion(1)/mudispersion(2);
				betaminusline = beta_arr - gd*omega_arr;
				[pdispersion2,~,mudispersion2]=polyfit(omega_arr,betaminusline,2);
				CD_mode_mm=-omega_arr(2)/(lambda_nm*1e-9)*2*pdispersion2(1)/(mudispersion2(2))^2;
				chromatic_dispersion_psnmkm(mm) = CD_mode_mm*scale_smm_to_psnmkm;
			end
			fib_params('CD_coeffs_psnmkm') = chromatic_dispersion_psnmkm;

			% store all parameters
			fib_params('D') = D;
			fib_params('fields') = fields; fib_params('fields_left') = fields_left; fib_params('fields_right') = fields_right;
			fib_params('neff') = neff; fib_params('neff_left') = neff_left; fib_params('neff_right') = neff_right;

			fib_params('MD_coeffs_psm') = modal_dispersion_coeffs_psm;
			fib_params('CD_coeffs_psnmkm') = chromatic_dispersion_psnmkm;
			
		end

		function c = get_speed_light()
			c = 299792458; % meters per second
		end

		function dn_out = apply_gaussian_filter_radial(dn_in, rho_arr, sigma_gaussian)
			% Apply a 1-D Gaussian filter with STD sigma_gaussian to the radial index distribution

			dr = rho_arr(2) - rho_arr(1);

			n_taps = ceil(6*sigma_gaussian/dr);
			if(mod(n_taps,2) == 0)
				n_taps = n_taps + 1; % need an odd number of taps
			end

			gaussian = @(sigma, x) (1/(sigma*sqrt(2*pi)))* exp(-0.5*x.^2 / sigma^2);
			x_taps = (-(n_taps-1)/2:(n_taps-1)/2)*dr;
			g_filt = gaussian(sigma_gaussian, x_taps);
			delay = (n_taps-1)/2;

			dn_out = filter(g_filt, 1, [dn_in zeros(1,delay)]) / sum(abs(g_filt));
			%post-process the delay
			dn_out = dn_out(delay+1:end);
		end

		function [dn_out, rho_out] = interpolate_index_radial(dn_in, rho_in, interp_factor)
			% Apply the interp1 function to test the effects of gridding
			% interp_factor = interpolation factor, < 1 if we want to shrink, > 1 if we want to expand
			n_pts_in = length(rho_in);
			n_pts_out = round(n_pts_in*interp_factor);
			rho_out = linspace(min(rho_in), max(rho_in), n_pts_out);
			dn_out = interp1(rho_in, dn_in, rho_out);
		end

		function print_fiber_summary()
			bank = fiber_bank();

			fib = bank.get_GI_MMF_1();
			fib = utils.solve_fiber_properties(fib);
			disp(sprintf('\nGI_MMF_1 @ %.2f nm\n', fib('center_wavelength_nm')));
			disp(sprintf('D = %d modes', fib('D')));
			disp(sprintf('rms group delay [ps/m] = %f', rms(fib('MD_coeffs_psm'))));
			disp(sprintf('rms chromatic dispersion [ps/nm*km] = %f', rms(fib('CD_coeffs_psnmkm'))));
			disp(sprintf('average modal effective area [um^2] = %f', mean(fib('Aeff'))));

			fib = bank.get_GI_MMF_2();
			fib = utils.solve_fiber_properties(fib);
			disp(sprintf('\nGI_MMF_2 @ %.2f nm\n', fib('center_wavelength_nm')));
			disp(sprintf('D = %d modes', fib('D')));
			disp(sprintf('rms group delay [ps/m] = %f', rms(fib('MD_coeffs_psm'))));
			disp(sprintf('rms chromatic dispersion [ps/nm*km] = %f', rms(fib('CD_coeffs_psnmkm'))));
			disp(sprintf('average modal effective area [um^2] = %f', mean(fib('Aeff'))));

			fib = bank.get_Corning_SMF28();
			fib = utils.solve_fiber_properties(fib);
			disp(sprintf('\nCorning SMF28 @ %.2f nm\n', fib('center_wavelength_nm')));
			disp(sprintf('D = %d modes', fib('D')));
			disp(sprintf('rms group delay [ps/m] = %f', rms(fib('MD_coeffs_psm'))));
			disp(sprintf('rms chromatic dispersion [ps/nm*km] = %f', rms(fib('CD_coeffs_psnmkm'))));
			disp(sprintf('average modal effective area [um^2] = %f', mean(fib('Aeff'))));

			fib = bank.get_StepIndex_MMF();
			fib = utils.solve_fiber_properties(fib);
			disp(sprintf('\nStepIndex_MMF @ %.2f nm\n', fib('center_wavelength_nm')));
			disp(sprintf('D = %d modes', fib('D')));
			disp(sprintf('rms group delay [ps/m] = %f', rms(fib('MD_coeffs_psm'))));
			disp(sprintf('rms chromatic dispersion [ps/nm*km] = %f', rms(fib('CD_coeffs_psnmkm'))));
			disp(sprintf('average modal effective area [um^2] = %f', mean(fib('Aeff'))));

			fib = bank.get_MMF_withTrench_1();
			fib = utils.solve_fiber_properties(fib);
			disp(sprintf('\nMMF with Trench JLT 2016 @ %.2f nm\n', fib('center_wavelength_nm')));
			disp(sprintf('D = %d modes', fib('D')));
			disp(sprintf('rms group delay [ps/m] = %f', rms(fib('MD_coeffs_psm'))));
			disp(sprintf('rms chromatic dispersion [ps/nm*km] = %f', rms(fib('CD_coeffs_psnmkm'))));
			disp(sprintf('average modal effective area [um^2] = %f', mean(fib('Aeff'))));
		end
		
	end
end