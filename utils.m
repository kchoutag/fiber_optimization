classdef utils
	methods(Static)
		function plot_cell_array(fig_name, cell_in, x_label, y_label)
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
			%dsfig(fig_name);
			plot(cell_array,'-b','linewidth',2);
			xlabel(x_label);
			ylabel(y_label);
		end

		function plot_results(fig_name, fiber_params)
			dsfig(fig_name);
			
			x_arr = linspace(-fiber_params('dx')*fiber_params('nx')*1e6/2, fiber_params('dx')*fiber_params('nx')*1e6/2, fiber_params('nx'));
			y_arr = linspace(-fiber_params('dy')*fiber_params('ny')*1e6/2, fiber_params('dy')*fiber_params('ny')*1e6/2, fiber_params('ny'));
			
			subplot(221);
			imagesc(x_arr, y_arr, fiber_params('index_distr')); colorbar;
            xlabel('x (\mu m)'); ylabel(' y (\mu m)');
		end
	end
end