function [ x ] = istft( stftx, stftshift, iscomplex )
% perform short-time Fourier transform
% stftx:            the stft data
% stftshift:        stft shift
% iscomplex:        true means the output is complex data
%
% x:            the transformed result
%

[F, T]=size(stftx);

% real or complex data
if iscomplex
    fftsize=F;
else 
    fftsize=(F-1)*2;
end

% synthesis window
% w=hann(fftsize);
w=hann(fftsize, 'periodic');
swin=sqrt(w*2.0*stftshift/fftsize);

% w=chebwin(fftsize);
% if length(w)==128
%     w=w/1.472905863110308;
% elseif length(w)==256
%     w=w/1.476632674000508;
% end
% 
% swin=sqrt(w);

% data length
len=dataLength(T, stftshift, fftsize);

%% perform istft
x=zeros(len, 1);
fx=zeros(fftsize, 1);
idx1=1;

for tau=1:T
    fx(1:F)=stftx(:, tau);
    % the complex conjugate part
    if ~iscomplex
        for fi=2:F-1
            fx(fftsize-fi+2)=conj(fx(fi));
        end
    end
    
    tx=ifft(fx, fftsize);
    tx=tx.*swin;
    if ~iscomplex
        tx=real(tx);
    end
    
    x(idx1:idx1+fftsize-1)=x(idx1:idx1+fftsize-1)+tx(:);
    
    idx1=idx1+stftshift;
end

end
