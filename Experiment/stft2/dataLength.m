function [ len ] = dataLength( tau, stftshift, fftsize )
% calculate the data length of the istft output
% tau:              no. of stft frames
% stftshift:        stft shift
% fftsize:          fft size
%
% len:              output data length
%

len=(tau-1)*stftshift+fftsize;

end
