clear; clc; close all;
tic;

gaussian = @(sigma, x) (1/(sigma*sqrt(2*pi)))* exp(-0.5*x.^2 / sigma^2);

sigma = 0.05;
x_filt = -0.75:0.05:0.75;
x_data = 0:0.15:100;

g = gaussian(sigma, x_filt);
delay = (length(x_filt)-1)/2;

%signal = randn(1,length(x_data));
signal = 0.2*sin(2*pi*0.1*x_data) + 0.05*randn(1,length(x_data));

% filter the signal
signal_filt = filter(g, 1, [signal zeros(1,delay-1)]) / sum(abs(g));

%post-process the delay
signal_filt = signal_filt(delay:end);

subplot(311);
plot(1:length(g),g); title('filter');
subplot(312);
plot(1:length(signal),signal); title('signal');
subplot(313);
plot(1:length(signal_filt), signal_filt); title('filtered');

toc;
