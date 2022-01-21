function waveform = sm_get_waveforms(rawtrace, spk)

fs = spk.fs;
spk = spk.spk;
nnum = length(spk);
waveform = cell(nnum,2);
for ii = 1:length(spk)
    spktimes = spk(ii).spiketimes;
    % convert spike times in ms to sampling index
    spktimes = round(spktimes * (fs/1000)); 
    % remove spikes occur within 2 ms towards the two ends of the recording
    spktimes(spktimes<=20) = [];
    spktimes(spktimes>=size(rawtrace,2)-20) = [];
    
    spkwaves = zeros(size(rawtrace, 1),41,length(spktimes));
    for kk = 1:length(spktimes)
        spkwaves(:,:,kk) = rawtrace(:,spktimes(kk)-20:spktimes(kk)+20);
    end
    
    spkwavavg = mean(spkwaves,3);
    spkwavstd = std(spkwaves,0, 3);
    waveform{ii,1} = spkwavavg;
    waveform{ii, 2} = spkwavstd;
end



end