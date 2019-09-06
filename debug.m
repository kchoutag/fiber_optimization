clear; clc;
bank = fiber_bank();
fiber_params = bank.get_GI_MMF_1();
fiber_params = bank.solve_fiber_properties(fiber_params);