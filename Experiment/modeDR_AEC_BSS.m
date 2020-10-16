function dataout=modeDR_AEC_BSS(nummics, numrefs, datain)
%
% Perform dr, aec, and bss.
% nummics:              no. of mic channels
% numrefs:              no. of reference channels
% datain:               input data
% dataout:              output data
%

%% perform stft
addpath('stft2');
% fft size
fftsize=512;
stftshift=fftsize/2;

M=nummics;
R=numrefs;
N=M;

Xtf=cell(M+R, 1);
for m=1:M+R
    Xtf{m}=stft(datain(:, m), stftshift, fftsize, false);
end
[K, T]=size(Xtf{1});

Ytf=cell(M, 1);
for m=1:M
    Ytf{m}=zeros(K, T);
end

%% params
% dr filter length
DR_FLEN=5;
% forgetting factor for dr and aec
DRAEC_FORGET=0.999;
% forgetting factor for bss
BF_FORGET=0.999;
% the shape parameter of the source prior
GAMMA=0.2;
%
% used to keep stable
%
VAR_BIAS=0.01;
STABLE_EPS=1e-3;
DRAEC_DIAGLOAD=1e-6;
BF_DIAGLOAD=1e-6;

%% space for dr
% current mic data backup
Miccurrent=zeros(K, M);
% delayed mic data
drfsize=M*DR_FLEN;
Micdelay=zeros(K, drfsize);

% mic and delay correlation
Cmd=cell(K, 1);
for k=1:K
    Cmd{k}=zeros(M, drfsize);
end

% delay auto correlation
Cdd=cell(K, 1);
for k=1:K
    Cdd{k}=zeros(drfsize, drfsize);
end

% reverb path
Reverbpath=cell(K, 1);
for k=1:K
    Reverbpath{k}=zeros(M, drfsize);
end

%% space for aec
% mic-reference correlation
Cmr=cell(K, 1);
for k=1:K
    Cmr{k}=zeros(M, R);
end

% reference auto correlation
Crr=cell(K, 1);
for k=1:K
    Crr{k}=zeros(R, R);
end

% echo path
Echopath=cell(K, 1);
for k=1:K
    Echopath{k}=zeros(M, R);
end

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

%% perform iteration
for tau=1:T    
    %% perform dr
    % direct and early reverberation
    Early=zeros(K, M);
    
    %
    % shift in new data
    %
    % shift old mic data back
    Micdelay=circshift(Micdelay, M, 2);
    % shift in delayed mic data
    Micdelay(:, 1:M)=Miccurrent;
    
    % mic data backup
    for m=1:M
        Miccurrent(:, m)=Xtf{m}(:, tau);
    end
    
    for k=1:K
        % calculate late reverberation
        ref=Micdelay(k, :).';
        late=Reverbpath{k}*ref;
        
        % direct and early reverberation
        mic=Miccurrent(k, :).';
        early=mic-late;
        % output data
        Early(k, :)=early.';
        
        %
        % calculate nonlinearity
        %
        xsq=abs(mic).^2;
        ysq=abs(early).^2;
        
        phi=0;
        for m=1:M
            if ysq(m)<=xsq(m)
                phi=phi+ysq(m);
            else
                phi=phi+xsq(m);
            end
        end
        
        phi=(1-DRAEC_FORGET)*(phi+VAR_BIAS)^((GAMMA-2)/2);
        
        % update mic ref correlation
        Cmd{k}=DRAEC_FORGET*Cmd{k}+phi*(mic*ref');
        
        % update ref auto-correlation
        Cdd{k}=DRAEC_FORGET*Cdd{k}+phi*(ref*ref');
        
        % update reverb path
        Reverbpath{k}=Cmd{k}/(Cdd{k}+DRAEC_DIAGLOAD*eye(drfsize, drfsize));
    end
    
    %% perform aec
    % nearend data
    Nearend=zeros(K, M);
    
    for k=1:K
        % reference data
        ref=zeros(R, 1);
        for r=1:R
            ref(r)=Xtf{M+r}(k, tau);
        end
        
        % echo
        echo=Echopath{k}*ref;
        
        % nearend
        mic=Early(k, :).';
        nearend=mic-echo;
        % output data
        Nearend(k, :)=nearend.';
        
        %
        % calculate nonlinearity
        %
        xsq=abs(mic).^2;
        ysq=abs(nearend).^2;
        
        phi=0;
        for m=1:M
            if ysq(m)<=xsq(m)
                phi=phi+ysq(m);
            else
                phi=phi+xsq(m);
            end
        end
        
        phi=(1-DRAEC_FORGET)*(phi+VAR_BIAS)^((GAMMA-2)/2);
        
        % update mic ref correlation
        Cmr{k}=DRAEC_FORGET*Cmr{k}+phi*(mic*ref');
        
        % update ref auto-correlation
        Crr{k}=DRAEC_FORGET*Crr{k}+phi*(ref*ref');
        
        % update echo path
        Echopath{k}=Cmr{k}/(Crr{k}+DRAEC_DIAGLOAD*eye(R, R));
    end
    
    %% perform bss
    Bssout=zeros(K, M);
    
    %
    % calculate nonlinearity
    %
    phi1=0;
    phi2=0;
    
    for k=1:K
        x=Nearend(k, :).';
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
        x=Nearend(k, :).';
        C1{k}=BF_FORGET*C1{k}+phi1*(x*x');
        C2{k}=BF_FORGET*C2{k}+phi2*(x*x');
        
        % solve gev problem
        D=heig2(BF_DIAGLOAD, C2{k}, C1{k});
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
