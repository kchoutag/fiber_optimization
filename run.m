clear; clc; close all;
tic;
if (feature('ShowFigureWindows') && strcmp(get(0,'DefaultFigureVisible'), 'on'))
    set(0,'DefaultFigureWindowStyle','docked');
end
set(0,'defaultaxesfontsize',14);
set(gca, 'FontName', 'Arial'); close(gcf);

test = experiment();
%test.test_robustness_GI_MMF_low_MD();
%results = test.MingJunLi_vs_Karthik_1(); % Ming Jun Li's fiber and his optimization
%results = test.MingJunLi_vs_Karthik_2();  % Karthik's optimization on Ming Jun Li's fiber
test.MMF_GI_radial_reduce_modal_dispersion(); % WORKING
%test.MMF_GI_radial_reduce_modal_dispersion_part2(); % 
%test.RCF_freeform_increase_degeneracies(); % WORKING
%test.MMF_SI_radial_reduce_coupling();
%test.MMF_GI_radial_reduce_coupling();
%test.SMF_SI_increase_CD();
%test.run_playground();


%utils.print_fiber_summary();


%bank = fiber_bank();
%foo = bank.get_Corning_SMF28();
%foo = utils.solve_fiber_properties(foo);
%utils.visualize_fiber(foo);

%lambda_arr = linspace(1500,1600,50)*1e-9;
%D_variation = test_robustness.vary_wavelength_number_modes(foo, lambda_arr);
%figure();
%plot(lambda_arr*1e9, D_variation, 'linewidth', 2);
%xlabel('Wavelength (nm)'); ylabel('Number of modes');


toc;