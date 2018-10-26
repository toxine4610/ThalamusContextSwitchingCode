function jittered_spikes = jitterspike(spiketimes,jitter);

x = rand(1,length(spiketimes)).*(jitter);
y = x - mean(x);
jittered_spikes = spiketimes + y';
