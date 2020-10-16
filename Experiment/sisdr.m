function val=sisdr(s, y)
%
% Calculate the scale-invariant SDR. 
%
% Le Roux, Jonathan, et al. "SDR¨Chalf-baked or well done?." ICASSP
% 2019-2019 IEEE International Conference on Acoustics, Speech and Signal
% Processing (ICASSP). IEEE, 2019.
% 
% s:                    source signal 
% y:                    estimated signals as column vectors
% val:                  the sisdr index
%

vals=zeros(size(y, 2), 1);
for n=1:length(vals)
    vals(n)=sisdr1(s, y(:, n));
end

val=max(vals);

end

function val=sisdr1(s, y)

ss=s*(y'*s/(s'*s));
res=ss-y;
val=10*log10(ss'*ss/(res'*res));

end
