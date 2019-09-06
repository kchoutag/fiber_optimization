classdef optimizers
	properties
	end
	methods
		function obj = optimizers()
		end

		function [d_index_distr] = opt_mode_coupling_freeform(obj, fiber_params, init_fiber_params, opt_params)
			D = fiber_params('D');
			nx = fiber_params('nx');
			ny = fiber_params('ny');

			init_neff = init_fiber_params('neff');
			curr_neff = fiber_params('neff');
			curr_fields = fiber_params('fields');

			% set the optimization target
			if(opt_params('direction') == 'MIN')
				des_neff = linspace(init_neff(1), init_neff(end), init_fiber_params('D'));
			elseif (opt_params('direction') == 'MAX')
				des_neff = ones(1,init_fiber_params('D'))*mean(init_neff);
            end

			d_index_distr = zeros(ny,nx);
			for xx = 1:nx
				for yy = 1:ny
					for mm = 1:min(init_fiber_params('D'), D)
						d_index_distr(yy,xx) = d_index_distr(yy,xx) - ((curr_neff(mm) - des_neff(mm))/(curr_neff(mm)))*curr_fields(yy,xx,mm)^2;
					end
				end
			end

			% normalize index update
			max_d_index = max(max(d_index_distr));
			d_index_distr = opt_params('max_dn')*d_index_distr/max_d_index;
		end


	end
end