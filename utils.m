classdef utils
	methods(Static)
		function plot_cell_array(fig_name, x_vals, cell_in, x_label, y_label)
			dsfig(fig_name);

			dim1 = length(cell_in);
			dim2 = 0;
			for ii = 1:dim1
				dim2 = max(dim2, length(cell_in{ii}));
			end
			cell_array = nan(dim1,dim2);
			for ii = 1:dim1
				cell_array(ii,1:length(cell_in{ii})) = cell_in{ii};
            end
            
            if(length(x_vals) > 1)
                plot(x_vals, cell_array,'-b','linewidth',2);
                xlabel(x_label);
                ylabel(y_label);
                axis tight;
            end
		end

		function plot_results(fig_name, fiber_params)
			dsfig(fig_name);
			
			x_arr = linspace(-fiber_params('dx')*fiber_params('nx')*1e6/2, fiber_params('dx')*fiber_params('nx')*1e6/2, fiber_params('nx'));
			y_arr = linspace(-fiber_params('dy')*fiber_params('ny')*1e6/2, fiber_params('dy')*fiber_params('ny')*1e6/2, fiber_params('ny'));
			
			subplot(221);
			n_clad = utils.get_index_at_wavelength(fiber_params('center_wavelength_nm'));
			n_xy = n_clad + fiber_params('index_distr_offset_from_cladding');
			imagesc(x_arr, y_arr, n_xy); axis square; colorbar;
            xlabel('x (\mum)'); ylabel(' y (\mum)');

            subplot(222);
            [X, Y] = meshgrid(x_arr, y_arr);
			surf(X, Y, n_xy); colorbar;
			axis square; shading interp;
			xlabel('x (\mum)');
			ylabel('y (\mum)');

			subplot(223);
			plot(fiber_params('Aeff'),'--o');
			xlabel('Mode index'); ylabel('Effective area (\mum^2)');
			axis tight;
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

		function [fib_params] = solve_fiber_properties(fib_params)
			addpath('../modesolver-2011-04-22/');
			n_modes_upper_lim = 70;
			c = utils.get_speed_light();
			lambda_nm = fib_params('center_wavelength_nm');
			d_lambda_nm = fib_params('d_wavelength_nm');

			% center wavelength
			n_clad = utils.get_index_at_wavelength(lambda_nm);
			n_xy = fib_params('index_distr_offset_from_cladding') + n_clad;
			n_core = max(max(n_xy)); 
			eps_profile = (n_xy).^2;
			[fields, neff_aux] = svmodes(lambda_nm*1e-9, n_core, n_modes_upper_lim, fib_params('dx'), fib_params('dy'), eps_profile, '0000', 'scalar');
			valid_idx = neff_aux > n_clad;
			fields = fields(:,:,valid_idx);
			neff = neff_aux(valid_idx);
			D_center = length(neff);

			% left wavelength
			left_wavelength_nm = lambda_nm - d_lambda_nm;
			n_clad = utils.get_index_at_wavelength(left_wavelength_nm);
			n_xy = fib_params('index_distr_offset_from_cladding') + n_clad;
			n_core = max(max(n_xy)); 
			eps_profile = (n_xy).^2;
			[fields_left, neff_aux_left] = svmodes(left_wavelength_nm*1e-9, n_core, n_modes_upper_lim, fib_params('dx'), fib_params('dy'), eps_profile, '0000', 'scalar');
			valid_idx = neff_aux_left > n_clad;
			fields_left = fields_left(:,:,valid_idx);
			neff_left = neff_aux_left(valid_idx);
			D_left = length(neff_left);

			% right wavelength
			right_wavelength_nm = lambda_nm + d_lambda_nm;
			n_clad = utils.get_index_at_wavelength(right_wavelength_nm);
			n_xy = fib_params('index_distr_offset_from_cladding') + n_clad;
			n_core = max(max(n_xy)); 
			eps_profile = (n_xy).^2;
			[fields_right, neff_aux_right] = svmodes(right_wavelength_nm*1e-9, n_core, n_modes_upper_lim, fib_params('dx'), fib_params('dy'), eps_profile, '0000', 'scalar');
			valid_idx = neff_aux_right > n_clad;
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

			% calculate effective areas
			Aeff = zeros(1,D);
			dx = fib_params('dx'); dy = fib_params('dy');
			for ii = 1:D
				% Aeff in terms of um^2
				Aeff(ii) = 1e6*1e6 * sum(sum(abs(fields(:,:,ii)).^2 * dx * dy))^2 / sum(sum(abs(fields(:,:,ii)).^4 * dx * dy));
			end
			fib_params('Aeff') = Aeff;

			% calculate the modal dispersion for all modes using the central difference approximation (in arbitrary units, TODO: fix units?)
    		modal_dispersion_coeffs = (neff_right - neff_left)/(2*fib_params('d_wavelength_nm')*1e-9);

    		% calculate chromatic dispersion for all modes (in arbitrary units)
    		%chromatic_dispersion_coeffs = (-fib_params('center_wavelength_nm')*1e-9/(c * (fib_params('d_wavelength_nm')*1e-9)^2)) * ...
			%   				  (neff_left + neff_right - 2*neff);

			k0 = 2*pi/(lambda_nm*1e-9);
			beta = k0*neff; beta_left = k0*neff_left; beta_right = k0*neff_right;
			d_omega = (-2*pi*c/((lambda_nm*1e-9)^2)) * d_lambda_nm * 1e-9;
			beta2 = (beta_left + beta_right - 2*beta)/(d_omega^2);
			chromatic_dispersion_coeffs = (-2*pi*c/((lambda_nm*1e-9)^2)) * beta2; % units ps/(nm*km) % TODO check units-how to convert to ps, nm and km?

			% store all parameters
			fib_params('D') = D;
			fib_params('fields') = fields; fib_params('fields_left') = fields_left; fib_params('fields_right') = fields_right;
			fib_params('neff') = neff; fib_params('neff_left') = neff_left; fib_params('neff_right') = neff_right;

			fib_params('modal_dispersion_coeffs') = modal_dispersion_coeffs;
			fib_params('chromatic_dispersion_coeffs') = chromatic_dispersion_coeffs;
			
		end

		function c = get_speed_light()
			c = 299792458; % meters per second
		end
		
	end
end