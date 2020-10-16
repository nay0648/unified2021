function dataout=modeBSS(nummics, numrefs, datain)
%
% Use Aux-IVA to solve dr, aec and separation.
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
% data buffer
fsize=M+R+M*DR_FLEN;
Bufdata=zeros(K, fsize);
% current mic backup
Miccurrent=zeros(K, M);

% the weighted correlation matrices
C1=cell(K, 1);
C2=cell(K, 1);
for k=1:K
    C1{k}=STABLE_EPS*eye(fsize, fsize);
    C2{k}=STABLE_EPS*eye(fsize, fsize);
end

% demixing matrices
Demix=cell(K, 1);
for k=1:K
    Demix{k}=eye(fsize, fsize);
end

%% perform iteration
for tau=1:T
    Bssout=zeros(K, M);
    
    %
    % shift in new data
    %
    % shift old mic data back
    Bufdata=circshift(Bufdata, M, 2);
    % shift in delayed mic data
    Bufdata(:, M+R+1:M+R+M)=Miccurrent;
    
    % mic data backup
    for m=1:M
        Miccurrent(:, m)=Xtf{m}(:, tau);
    end
    % current mic data
    Bufdata(:, 1:M)=Miccurrent;
    
    % shift in reference data
    for r=1:R
        Bufdata(:, M+r)=Xtf{M+r}(:, tau);
    end
    
    %
    % calculate nonlinearity
    %
    phi1=0;
    phi2=0;
    
    for k=1:K
        x=Bufdata(k, :).';
        y=Demix{k}*x;
        % output data
        Bssout(k, :)=y(1:M).';
        
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
        x=Bufdata(k, :).';
        C1{k}=BF_FORGET*C1{k}+phi1*(x*x');
        C2{k}=BF_FORGET*C2{k}+phi2*(x*x');
        
        %
        % update demixing matrix
        %
        H1=inv(Demix{k}*C1{k}+BF_DIAGLOAD*eye(fsize, fsize));
        w1=H1(:, 1);
        
        H2=inv(Demix{k}*C2{k}+BF_DIAGLOAD*eye(fsize, fsize));
        w2=H2(:, 2);
        
        D=eye(fsize, fsize);
        D(1, :)=w1';
        D(2, :)=w2';
        
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
        
        D=eye(fsize, fsize);
        D(1, :)=a1*w1';
        D(2, :)=a2*w2';
        
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
