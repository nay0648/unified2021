function [noise2, gain] = adjustSNR(signal, noise, snrdb)
%
% Adjust noise signal amplitude to meet the specified SNR.
% signal:               target signals
% noise:                noise signals
% snrdb:                required SNR (dB)
% noise2:               adjusted noise signals
%

blocksize=320;
% background threshold
noisefloordb=-55;

%% calculate signal power
totalbpow=0;
numvalidblocks=0;
for bi=1:floor(size(signal, 1)/blocksize)
    block=signal((bi-1)*blocksize+1:(bi-1)*blocksize+blocksize, :);
    bsq=abs(block).^2;
    bpow=sum(bsq(:))/(size(bsq, 1)*size(bsq, 2));
    
    bpowdb=10*log10(bpow);
    if bpowdb>=noisefloordb
        totalbpow=totalbpow+bpow;
        numvalidblocks=numvalidblocks+1;
    end
end

pows=totalbpow/numvalidblocks;

%% calculate noise power
totalbpow=0;
numvalidblocks=0;
for bi=1:floor(size(noise, 1)/blocksize)
    block=noise((bi-1)*blocksize+1:(bi-1)*blocksize+blocksize, :);
    bsq=abs(block).^2;
    bpow=sum(bsq(:))/(size(bsq, 1)*size(bsq, 2));
    
    bpowdb=10*log10(bpow);
    if bpowdb>=noisefloordb
        totalbpow=totalbpow+bpow;
        numvalidblocks=numvalidblocks+1;
    end
end

pown=totalbpow/numvalidblocks;

%% adjust snr
snr=10^(snrdb/10);
gain=sqrt(pows/(pown*snr));
noise2=gain*noise;

end
