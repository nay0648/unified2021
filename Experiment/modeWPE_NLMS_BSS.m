function dataout=modeWPE_NLMS_BSS(nummics, numrefs, datain)
%
% Perform wpe, nlms, and bss.
% nummics:              no. of mic channels
% numrefs:              no. of reference channels
% datain:               input data
% dataout:              output data
%

fs=16000;

%% perform wpe
addpath('wpe_v1.33');

mic=datain(:, 1:2);
ref=datain(:, 3);

cfgs='wpe_v1.33/settings/local.m';
early=wpe(mic, cfgs); 

%% perform nlms
addpath('Speex-AEC-matlab-master');

%    Usage: 
%
%       speex_mdf_out = speex_mdf(Fs, u, d, filter_length, frame_size, dbg_var_name);
%       
%       Fs                  sample rate
%       u                   speaker signal, column vector in range [-1; 1]
%       d                   microphone signal, column vector in range [-1; 1]
%       filter_length       typically 250ms, i.e. 4096 @ 16k FS 
%                           must be a power of 2
%       frame_size          typically 8ms, i.e. 128 @ 16k Fs 
%                           must be a power of 2
%       dbg_var_name        internal state variable name to trace. 
%                           Default: 'st.leak_estimate'.
%
%    Jonathan Rouach <jonr@waves.com>
%    
filter_length=512;
frame_size=256;
nearend1=speex_mdf(fs, ref, early(:, 1), filter_length, frame_size);
nearend2=speex_mdf(fs, ref, early(:, 2), filter_length, frame_size);
nearend=[nearend1.e, nearend2.e];

%% perform stft
addpath('stft2');
% fft size
fftsize=512;
stftshift=fftsize/2;

M=nummics;
N=M;

Xtf=cell(M, 1);
for m=1:M
    Xtf{m}=stft(nearend(:, m), stftshift, fftsize, false);
end
[K, T]=size(Xtf{1});

Ytf=cell(M, 1);
for m=1:M
    Ytf{m}=zeros(K, T);
end

%% params
% forgetting factor for bss
BF_FORGET=0.999;
% the shape parameter of the source prior
GAMMA=0.2;
%
% used to keep stable
%
VAR_BIAS=0.01;
STABLE_EPS=1e-3;
BF_DIAGLOAD=1e-6;

%% space for bss
% the weighted correlation matrices
C1=cell(K, 1);
C2=cell(K, 1);
for k=1:K
    C1{k}=STABLE_EPS*eye(M, M);
    C2{k}=STABLE_EPS*eye(M, M);
end

% demixing matrices
Demix=cell(K, 1);
for k=1:K
    Demix{k}=eye(N, M);
end

%% perform bss
for tau=1:T
    %% perform bss
    Bssout=zeros(K, M);
    
    %
    % calculate nonlinearity
    %
    phi1=0;
    phi2=0;
    
    for k=1:K
        x=zeros(M, 1);
        for m=1:M
            x(m)=Xtf{m}(k, tau);
        end
        
        y=Demix{k}*x;
        % output data
        Bssout(k, :)=y.';
        
        phi1=phi1+abs(y(1))^2;
        phi2=phi2+abs(y(2))^2;
    end
    
    phi1=(1-BF_FORGET)*(phi1+VAR_BIAS)^((GAMMA-2)/2);
    phi2=(1-BF_FORGET)*(phi2+VAR_BIAS)^((GAMMA-2)/2);
    
    % update the demixing matrices
    for k=1:K
        %
        % accumulate the weighted correlation
        %
        x=zeros(M, 1);
        for m=1:M
            x(m)=Xtf{m}(k, tau);
        end
        
        C1{k}=BF_FORGET*C1{k}+phi1*(x*x');
        C2{k}=BF_FORGET*C2{k}+phi2*(x*x');
        
        %
        % solve gev problem
        %
        [Ev, Ed]=eig(C2{k}+BF_DIAGLOAD*eye(M, M), C1{k}+BF_DIAGLOAD*eye(M, M));
        if Ed(1, 1)>=Ed(2, 2)
            e1=Ev(:, 1);
            e2=Ev(:, 2);
        else 
            e1=Ev(:, 2);
            e2=Ev(:, 1);
        end
    
        D=[e1'; e2'];
        
        %
        % solve the scaling ambiguity
        %
        A=inv(D);
        
        if abs(A(1, 1))>=abs(A(2, 1))
            a1=A(1, 1);
        else
            a1=A(2, 1);
        end
        
        if abs(A(2, 2))>=abs(A(1, 2))
            a2=A(2, 2);
        else
            a2=A(1, 2);
        end
        
        D=diag([a1; a2])*D;
        Demix{k}=D;
    end
    
    for m=1:M
        Ytf{m}(:, tau)=Bssout(:, m);
    end
end

%% perform istft and output signal
dataout=zeros(dataLength(T, stftshift, fftsize ), N);
for n=1:N
    dataout(:, n)=istft(Ytf{n}, stftshift, false);
end

end
