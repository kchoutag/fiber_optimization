clear; clc; close all;
tic;
if (feature('ShowFigureWindows') && strcmp(get(0,'DefaultFigureVisible'), 'on'))
    set(0,'DefaultFigureWindowStyle','docked');
end

test = experiment();
test.MMF_GI_radial_reduce_modal_dispersion();
%test.MMF_GI_freeform_reduce_coupling();
%test.run_playground();

%{
utils.print_fiber_summary();

bank = fiber_bank();
foo = bank.get_MMF_withTrench_1();
foo = utils.solve_fiber_properties(foo);
utils.visualize_fiber(foo);
%}


toc;