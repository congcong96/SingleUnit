function spk = calc_ISI_violation_percentage(spk, isivpthresh)

% calculate percentage of ISIs that violate the refractory period of the
% neuron. Default threshold is 2ms. Recommended 1.5ms?

if ~exist('isivp','var')
    isivpthresh = 2;
end

for i = 1:length(spk)
    try
        spktimes = spk(i).spiketimes;
    catch
        spktimes = spk(i).spiketimes_dmr;
    end
    diffst = diff(spktimes);
    spk(i).ISI_vio = sum(diffst <= isivpthresh) ./ length(diffst) * 100;
end

