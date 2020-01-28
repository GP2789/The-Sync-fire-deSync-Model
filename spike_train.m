function spikes = spike_train( n_neurons ,spikesPerS, duration_s, d_T)
% takes a random selection of spike trains from initial trains. Assumes
% that neurons in this group recieve noise from specific neurons.

n_trains = 10;
spikes = zeros(n_neurons, duration_s);
n_noise = n_trains * 10;
vt = (spikesPerS*(0.001*d_T)/n_trains) > rand(n_noise, duration_s); % make random numbers
rand_n = ceil(rand(n_trains,1,n_neurons)*size(vt,1));

for n=1:n_neurons
    spikes(n,:) = sum(vt(rand_n(:,1,n),:),1);
end

end