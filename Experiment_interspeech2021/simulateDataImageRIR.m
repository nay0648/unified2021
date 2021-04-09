function [testdata, target_interf_echo] = simulateDataImageRIR(sirdb, serdb, rt60, repeat, seed)
%
% Generate test data with simulated rir.
% sirdb:                specified SIR (dB)
% serdb:                specified SER (dB)
% rt60:                 reverberation time (s)
% testdata:             the generated data
% target_interf_echo:   target, interference, and echo signals
%

addpath('RIR-Generator-master');
% to reproduce experiments
rng(seed);

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
siglen=length(target);
% set target signal activity
target=[zeros(siglen, 1);target;target;zeros(siglen, 1)];

%
% interference
%
sub=dir(base_interf);
idx=randperm(size(sub, 1), 1);
interf=audioread([sub(idx).folder, '/', sub(idx).name]);
% set interference activity
begin=randi(length(interf)-4*siglen);
interf=interf(begin+1:begin+4*siglen);

%
% reference
%
sub=dir(base_ref);
idx=randperm(size(sub, 1), 1);
ref=audioread([sub(idx).folder, '/', sub(idx).name]);
% set reference activity
begin=randi(length(ref)-2*siglen);
ref=ref(begin+1:begin+2*siglen);
ref=[zeros(length(ref), 1);ref];

% repeat the same utterance
if repeat > 1
    target_bk=target; interf_bk=interf; ref_bk=ref;
    for i=1:repeat-1
        target = [target;target_bk];
        interf = [interf;interf_bk];
        ref = [ref;ref_bk];
    end
end

%% generate simulated rirs
room_bound_x=[3, 6];
room_bound_y=[4, 8];
room_bound_z=[2.5, 4];
% sample the room size
room_x=round(rand*(room_bound_x(2)-room_bound_x(1))+room_bound_x(1), 2);
room_y=round(rand*(room_bound_y(2)-room_bound_y(1))+room_bound_y(1), 2);
room_z=round(rand*(room_bound_z(2)-room_bound_z(1))+room_bound_z(1), 2);

c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
L = [room_x, room_y, room_z];  % Room dimensions [x y z] (m)
beta = rt60;                % Reverberation time (s)
n = 4096;                   % Number of samples
mtype = 'omnidirectional';  % Type of microphone
order = -1;                 % -1 equals maximum reflection order!
dim = 3;                    % Room dimension
orientation = 0;            % Microphone orientation (rad)
hp_filter = 1;              % Enable high-pass filter
early_reverb = fs * 0.05;   % define early reverberation

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
dist_bound = [0.2, 3];  % allowed source-receiver distance
while 1
    elevation=asin(2*rand - 1);
    azimuth=2*pi*rand;
    % favor large distance
    radii=dist_bound(1)+(dist_bound(2)-dist_bound(1))*(rand.^(1/3));
    [offset_x, offset_y, offset_z] = sph2cart(azimuth, elevation, radii);
    offset_xyz = [offset_x offset_y offset_z];
    s1 = offset_xyz + r0;
    % check if the source position is within the correct range, otherwise resample
    if s1 <= (L-space) & s1 >= space
        break
    end
end
while 1
    elevation=asin(2*rand - 1);
    azimuth=2*pi*rand;
    % favor large distance
    radii=dist_bound(1)+(dist_bound(2)-dist_bound(1))*(rand.^(1/3));
    [offset_x, offset_y, offset_z] = sph2cart(azimuth, elevation, radii);
    offset_xyz = [offset_x offset_y offset_z];
    s2 = offset_xyz + r0;
    % check if the source position is within the correct range, otherwise resample
    if s2 <= (L-space) & s2 >= space
        break
    end
end

% % show mic-source positions
% plot(r(:,1),r(:,2),'o'); xlim([0,L(1)]); ylim([0,L(2)]);
% hold on; plot(s1(1),s1(2), 'x'); plot(s2(1),s2(2), 'v'); hold off;

% loudspeaker is 15 cm below mic
s3=r0;
s3(3)=s3(3)-0.15;

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
target_early_reverb=conv(target, a11(1:early_reverb));
target_early_reverb(length(target)+1:end, :)=[];

interf_img=[conv(interf, a12), conv(interf, a22)];
interf_img(length(interf)+1:end, :)=[];
interf_early_reverb=conv(interf, a12(1:early_reverb));
interf_early_reverb(length(interf)+1:end, :)=[];

echo=[conv(ref, b1), conv(ref, b2)];
echo(length(ref)+1:end, :)=[];

%
% adjust snr
%
[target_img, gain]=adjustSNR(target, target_img, 0);
target_early_reverb=target_early_reverb*gain;
[interf_img, gain]=adjustSNR(target_img, interf_img, sirdb);
interf_early_reverb=interf_early_reverb*gain;
[echo, gain]=adjustSNR(target_img, echo, serdb);

% generate the mixed signals
testdata=[target_img+interf_img+echo, ref];
target_interf_echo=[target_early_reverb, interf_early_reverb, echo(:, 1)];

end
