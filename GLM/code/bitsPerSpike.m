function bps = bitsPerSpike(lambda, r, rbar)
% Bits per spike calculation
% bps = bitsPerSpike(lambda, r, rbar)
% must be in count per bin, not per second units
if numel(rbar)==1
    rbar=rbar*ones(size(lambda));
end

llsp   = r'*log(lambda) - sum(lambda);
llbase = r'*log(rbar) - sum(rbar);

bps = (llsp-llbase)/sum(r)/log(2);