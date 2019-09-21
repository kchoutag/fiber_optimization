clear; clc; close all;
tic;
if (feature('ShowFigureWindows') && strcmp(get(0,'DefaultFigureVisible'), 'on'))
    set(0,'DefaultFigureWindowStyle','docked');
end
set(0,'defaultaxesfontsize',14);
set(gca, 'FontName', 'Arial');

test = experiment();
test.MMF_GI_radial_reduce_modal_dispersion();
%test.MMF_SI_radial_reduce_coupling();
%test.SMF_SI_increase_CD();
%test.run_playground();

%{
utils.print_fiber_summary();

bank = fiber_bank();
foo = bank.get_MMF_withTrench_1();
foo = utils.solve_fiber_properties(foo);
lambda_arr = linspace(1500,1600,50)*1e-9;
D_variation = test_robustness.vary_wavelength_number_modes(foo, lambda_arr);
figure();
plot(lambda_arr*1e9, D_variation, 'linewidth', 2);
xlabel('Wavelength (nm)'); ylabel('Number of modes');
%}

toc;