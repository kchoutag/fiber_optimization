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
			imagesc(x_arr, y_arr, fiber_params('index_distr')); colorbar;
            xlabel('x (\mum)'); ylabel(' y (\mum)');

            subplot(222);
            [X, Y] = meshgrid(x_arr, y_arr);
			surf(X, Y, fiber_params('index_distr')); colorbar;
			axis square; shading interp;
			xlabel('x (\mum)');
			ylabel('y (\mum)');

			subplot(223);
			plot(fiber_params('Aeff'),'--o');
			xlabel('Mode index'); ylabel('Effective area (\mum^2)');
			axis tight;
		end
	end
end