function [ stftx ] = stft( x, stftshift, fftsize, iscomplex )
% perform short-time Fourier transform
% x:                the data
% stftshift:        stft shift
% fftsize:          fft size
% iscomplex:        true means the input is complex data
%
% stftx:            the transformed result
%

% analysis window
% w=hann(fftsize);
w=hann(fftsize, 'periodic');
awin=sqrt(w*2.0*stftshift/fftsize);

% w=chebwin(fftsize);
% if length(w)==128
%     w=w/1.472905863110308;
% elseif length(w)==256
%     w=w/1.476632674000508;
% end

% awin=sqrt(w);

% real or complex data
if iscomplex
    F=fftsize;
else 
    F=fftsize/2+1;
end

% no. of frames
T=numFrames(length(x), stftshift);

%% perform stft
stftx=zeros(F, T)+zeros(F, T)*1i;
idx1=1;

for tau=1:T
    tx=zeros(fftsize, 1);
    idx2=min(idx1+fftsize-1, length(x));
    
    tx(1:idx2-idx1+1)=x(idx1:idx2);
    tx=tx.*awin;
    
    fx=fft(tx, fftsize);
    stftx(:, tau)=fx(1:F);
    
    idx1=idx1+stftshift;
end

end
