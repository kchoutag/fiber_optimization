clear; clc;
%bank = fiber_bank();
%fiber_params = bank.get_GI_MMF_1();
%fiber_params = bank.solve_fiber_properties(fiber_params);
c_speed_light = 299792458;

B1 = 0.6961663;
B2 = 0.4079426;
B3 = 0.8974794;

L1 = 0.0684043; %um
L2 = 0.1162414; %um
L3 = 9.896161;  %um

omega1 = 2*pi*c_speed_light/(L1*1e-6);
omega2 = 2*pi*c_speed_light/(L2*1e-6);
omega3 = 2*pi*c_speed_light/(L3*1e-6);

wavelength_nm = linspace(500, 1600, 100);
n_fusedSi = zeros(1, length(wavelength_nm));

for ii = 1:length(wavelength_nm)
    omega_ii = 2*pi*c_speed_light/(wavelength_nm(ii)*1e-9);
    
    n_fusedSi(ii) = sqrt(1 + B1*((omega1^2)/(omega1^2 - omega_ii^2)) + ...
        B2*((omega2^2)/(omega2^2 - omega_ii^2)) + ...
        B3*((omega3^2)/(omega3^2 - omega_ii^2))); 
end

plot(wavelength_nm, n_fusedSi); xlabel('Wavelength (nm)'); ylabel('n fused Si');