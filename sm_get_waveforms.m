function waveform = sm_get_waveforms(rawtrace, spk)

fs = spk.fs;
spk = spk.spk;
nnum = length(spk);
waveform = cell(nnum,3);
hw = 1.5;
taxis = -hw:(1/fs*1000):hw;
nhw = hw/(1/fs*1000);
for ii = 1:length(spk)
    spktimes = spk(ii).spiketimes;
    % convert spike times in ms to sampling index
    spktimes = round(spktimes * (fs/1000)); 
    % remove spikes occur within 1.5 ms towards the two ends of the recording
    spktimes(spktimes<=nhw) = [];
    spktimes(spktimes>=size(rawtrace,2)-nhw) = [];
    
    spkwaves = zeros(size(rawtrace, 1),length(taxis),length(spktimes));
    for kk = 1:length(spktimes)
        spkwaves(:,:,kk) = rawtrace(:,spktimes(kk)-nhw:spktimes(kk)+nhw);
    end
    
    spkwavavg = mean(spkwaves,3);
    spkwavstd = std(spkwaves,0, 3);
    waveform{ii,1} = spkwavavg;
    waveform{ii, 2} = spkwavstd;
    waveform{ii, 3} = taxis;
end



end