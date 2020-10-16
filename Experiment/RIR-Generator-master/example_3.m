clear;
addpath('../');

c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
r = [2.5-0.05 1 1 ; 2.5+0.05 1 1];    % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)
s = [1.5 1 1];              % Source position [x y z] (m)
L = [5 7 2.4];                % Room dimensions [x y z] (m)
beta = 0.32;                 % Reverberation time (s)
n = 8192;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                 % -1 equals maximum reflection order!
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

h = rir_generator(c, fs, r, s, L, beta, n, mtype, order, dim, orientation, hp_filter);
h=10*h;
t=0:size(h, 2)-1;

figure(1);
plot(t, h(1, :), t, h(2, :));

figure(2);
rt60(h(1, :), fs)
