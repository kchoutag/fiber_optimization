clear; clc; close all;
tic;
if (feature('ShowFigureWindows') && strcmp(get(0,'DefaultFigureVisible'), 'on'))
    set(0,'DefaultFigureWindowStyle','docked');
end

test = experiment();
test.run_MMF_GI_reduce_coupling();
%test.run_playground();

toc;