function [ tau ] = numFrames( len, stftshift )
% calculate the no. of stft frames according to the input data length
% len:              input data length
% stftshift:        stft shift
%
% tau:              no. of stft frames
%   

tau=floor(len/stftshift);

end
