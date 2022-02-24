function spk = calc_firing_rate_stability(spk, trigger)
fs = spk(1).fs;
start_time = trigger(1)/fs*1000;
end_time = trigger(end)/fs*1000;
half_time = (start_time + end_time)/2;

for ii = 1:length(spk)
    try
        spiketimes = spk(ii).spiketimes;
    catch
        spiketimes = spk(ii).spiketimes_dmr;
    end
    spknum1 = sum(spiketimes > start_time & spiketimes < half_time);
    spknum2 = sum(spiketimes > half_time & spiketimes < end_time);
    stability = 1-abs((spknum2-spknum1)/(spknum1+spknum2));
    spk(ii).frs = stability;
end