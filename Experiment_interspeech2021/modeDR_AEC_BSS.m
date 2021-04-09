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
config;

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

%% params go to config

%% space for dr
% current mic data backup
Micbuffer=zeros(K, M*(DR_DELAY+1));
% delayed mic data
drfsize=M*DR_FLEN;
BufDR=zeros(K, drfsize);

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
aecfsize=R*AEC_FLEN;
BufAEC=zeros(K, aecfsize);
% mic-reference correlation
Cmr=cell(K, 1);
for k=1:K
    Cmr{k}=zeros(M, aecfsize);
end

% reference auto correlation
Crr=cell(K, 1);
for k=1:K
    Crr{k}=zeros(R, aecfsize);
end

% echo path
Echopath=cell(K, 1);
for k=1:K
    Echopath{k}=zeros(M, aecfsize);
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
for t=1:T    
    %% perform dr
    % direct and early reverberation
    Early=zeros(K, M);
    
    %
    % shift in mic data
    %
    Micbuffer=circshift(Micbuffer, M, 2);
    for m=1:M
        Micbuffer(:, m)=Xtf{m}(:, t);
    end
    % shift old mic data back
    BufDR=circshift(BufDR, M, 2);
    % shift in delayed mic data
    BufDR(:, 1:M)=Micbuffer(:,end-M+1:end);
    
    for k=1:K
        mic=permute(Micbuffer(k, 1:M), [2,1]);
        ref=permute(BufDR(k, :), [2,1]);
        
        % calculate late reverberation
        late=Reverbpath{k}*ref;
        
        % direct and early reverberation
        early=mic-late;
        % output data
        Early(k, :)=early;
        
        %
        % calculate nonlinearity
        %
        xsq=abs(mic).^2;
        ysq=abs(early).^2;
        
        phi=sum(ysq(ysq<xsq)) + sum(xsq(ysq>=xsq));
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
    
    BufAEC=circshift(BufAEC, R, 2);
    for r=1:R
        BufAEC(:,r)=Xtf{M+r}(:, t);
    end
    for k=1:K   
        mic=permute(Early(k, :), [2,1]);
        ref=permute(BufAEC(k, :), [2,1]);
        
        % echo
        echo=Echopath{k}*ref;
        
        % nearend
        nearend=mic-echo;
        % output data
        Nearend(k, :)=nearend;
        
        %
        % calculate nonlinearity
        %
        xsq=abs(mic).^2;
        ysq=abs(nearend).^2;
        
        phi=sum(ysq(ysq<xsq)) + sum(xsq(ysq>=xsq));
        phi=(1-DRAEC_FORGET)*(phi+VAR_BIAS)^((GAMMA-2)/2);
        
        % update mic ref correlation
        Cmr{k}=DRAEC_FORGET*Cmr{k}+phi*(mic*ref');
        
        % update ref auto-correlation
        Crr{k}=DRAEC_FORGET*Crr{k}+phi*(ref*ref');
        
        % update echo path
        Echopath{k}=Cmr{k}/(Crr{k}+DRAEC_DIAGLOAD*eye(aecfsize, aecfsize));
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
        Ytf{m}(:, t)=Bssout(:, m);
    end
end

%% perform istft and output signal
dataout=zeros(dataLength(T, stftshift, fftsize ), N);
for n=1:N
    dataout(:, n)=istft(Ytf{n}, stftshift, false);
end

end
