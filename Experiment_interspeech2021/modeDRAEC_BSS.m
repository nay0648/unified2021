function dataout=modeDRAEC_BSS(nummics, numrefs, datain)
%
% Perform dr and aec together, then bss.
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

%% space for dr and aec
% current mic data backup
Micbuffer=zeros(K, M*(DR_DELAY+1));
% ref and delayed mic data
draecfsize=R*AEC_FLEN+M*DR_FLEN;
Refmicdelay=zeros(K, draecfsize);

% mic-ref correlation
Cxr=cell(K, 1);
for k=1:K
    Cxr{k}=zeros(M, draecfsize);
end

% reference auto correlation
Crr=cell(K, 1);
for k=1:K
    Crr{k}=zeros(draecfsize, draecfsize);
end

% reverb and echo path
REPath=cell(K, 1);
for k=1:K
    REPath{k}=zeros(M, draecfsize);
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
    %% perform dr and aec together
    % direct nearend and early reverberation
    Early=zeros(K, M);
    
    %
    % shift in new data
    %
    Micbuffer=circshift(Micbuffer, M, 2);
    for m=1:M
        Micbuffer(:, m)=Xtf{m}(:, t);
    end
    
    % shift in reference data
    Refmicdelay(:, 1:R*AEC_FLEN)=circshift(Refmicdelay(:, 1:R*AEC_FLEN), R, 2);
    for r=1:R
        Refmicdelay(:, r)=Xtf{M+r}(:, t);
    end
    
    % delayed mic data
    Refmicdelay(:, R*AEC_FLEN+1:end)=circshift(Refmicdelay(:, R*AEC_FLEN+1:end), M, 2);
    Refmicdelay(:, R*AEC_FLEN+1:R*AEC_FLEN+M)=Micbuffer(:, end-M+1:end);
    
    for k=1:K
        mic=permute(Micbuffer(k, 1:M), [2,1]);
        ref=permute(Refmicdelay(k, :), [2,1]);

        % calculate late reverberation and echo
        late=REPath{k}*ref;
        
        % direct nearend and early reverberation
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
        Cxr{k}=DRAEC_FORGET*Cxr{k}+phi*(mic*ref');
        
        % update ref auto-correlation
        Crr{k}=DRAEC_FORGET*Crr{k}+phi*(ref*ref');
        
        % update echo and reverb path
        REPath{k}=Cxr{k}/(Crr{k}+DRAEC_DIAGLOAD*eye(draecfsize, draecfsize));
    end
    
    %% perform bss
    Bssout=zeros(K, M);
    
    %
    % calculate nonlinearity
    %
    phi1=0;
    phi2=0;
    
    for k=1:K
        x=Early(k, :).';
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
        x=Early(k, :).';
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
