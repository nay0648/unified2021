function [testdata, target_interf_echo] = simulateDataImageRIR(snrdb, serdb, rt60)
%
% Generate test data with simulated rir.
% snrdb:                specified SNR (dB)
% serdb:                specified SER (dB)
% rt60:                 reverberation time (s)
% testdata:             the generated data
% target_interf_echo:   target, interference, and echo signals
%

addpath('RIR-Generator-master');

%% base directories
base_target='target/*.wav';
base_interf='interference/*.wav';
base_ref='reference/*.wav';

%% load sources
%
% target
%
sub=dir(base_target);
idx=randperm(size(sub, 1), 1);
target=audioread([sub(idx).folder, '/', sub(idx).name]);

%
% interference
%
sub=dir(base_interf);
idx=randperm(size(sub, 1), 1);
interf=audioread([sub(idx).folder, '/', sub(idx).name]);

%
% reference
%
sub=dir(base_ref);
idx=randperm(size(sub, 1), 1);
ref=audioread([sub(idx).folder, '/', sub(idx).name]);

%% generate simulated rirs
c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
L = [5, 7, 2.4];            % Room dimensions [x y z] (m)
beta = rt60;                % Reverberation time (s)
n = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                 % -1 equals maximum reflection order!
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter

%
% generate mic position
%
d = 0.1;                    % mic spacing (m)
theta=rand(1, 1)*2*pi;
R=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
m1=R*[-d/2; 0; 0];
m2=R*[d/2; 0; 0];

space=0.5;
r0=space*ones(1, 3)+rand(1, 3).*(L-2*space*ones(1, 3));
r = [r0+m1'; r0+m2'];       % Receiver positions [x_1 y_1 z_1 ; x_2 y_2 z_2] (m)

%
% generate source position
%
space=0.1;
s1=space*ones(1, 3)+rand(1, 3).*(L-2*space*ones(1, 3));
space=0.2;
s2=space*ones(1, 3)+rand(1, 3).*(L-2*space*ones(1, 3));

% loudspeaker is 10 cm below mic
s3=r0;
s3(3)=s3(3)-0.1;

%
% generate rir
%
h=rir_generator(c, fs, r, s1, L, beta, n, mtype, order, dim, orientation, hp_filter);
a11=h(1, :);
a21=h(2, :);
% tt=0:length(a11)-1;
% plot(tt, a11, tt, a21);

h=rir_generator(c, fs, r, s2, L, beta, n, mtype, order, dim, orientation, hp_filter);
a12=h(1, :);
a22=h(2, :);

h=rir_generator(c, fs, r, s3, L, beta, n, mtype, order, dim, orientation, hp_filter);
b1=h(1, :);
b2=h(2, :);
% tt=0:length(b1)-1;
% plot(tt, b1, tt, b2);

%% mix data
%
% generate data image
%
target_img=[conv(target, a11), conv(target, a21)];
target_img(length(target)+1:end, :)=[];

interf_img=[conv(interf, a12), conv(interf, a22)];
interf_img(length(interf)+1:end, :)=[];

echo=[conv(ref, b1), conv(ref, b2)];
echo(length(ref)+1:end, :)=[];

%
% adjust snr
%
target_img=adjustSNR(target, target_img, 0);
interf_img=adjustSNR(target_img, interf_img, snrdb);
echo=adjustSNR(target_img, echo, serdb);

% generate the mixed signals
testdata=[target_img+interf_img+echo, ref];
target_interf_echo=[target_img(:, 1), interf_img(:, 1), echo(:, 1)];

end
