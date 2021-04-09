fs=16000;

% fft
fftsize = 1024;
stftshift = fftsize/2;

% aec filter length
AEC_FLEN = 5;
% dr filter delay
DR_DELAY = 2;
% dr filter length
DR_FLEN = 5;
% forgetting factor for dr and aec
DRAEC_FORGET=0.9995;
% forgetting factor for bss
BF_FORGET = 0.999;
% the shape parameter of the source prior
GAMMA = 0.2;
%
% used to keep stable
%
VAR_BIAS=0.01;
STABLE_EPS=1e-3;
DRAEC_DIAGLOAD=1e-6;
BF_DIAGLOAD=1e-6;
