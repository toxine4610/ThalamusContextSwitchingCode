function llsp = logliPoisson(lambda, r)
% Poisson log likelihood
% bps = bitsPerSpike(lambda, r)
% must be in count per bin, not per second units

etol = 1e-100;
lambda(lambda<etol)=etol;

llsp   = r'*log(lambda) - sum(lambda);