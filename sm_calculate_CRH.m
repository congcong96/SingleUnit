function [mtfhist, tmfaxis, smfaxis] = sm_calculate_CRH(spktrain, mtffile)
% calculate CRH with stim_mst generated mtf stim file and 0.5ms binned
% spike trains
% adapted from ne_calc_exp_site_nedata.mtfhist
% Input:
%   spktrain: nxm matrix, n neurons, m time bins
%   mtffile: mtf of stimulus (use batch_stimulus_to_tmf_smf to generate)
% Output:
% mtfhist: CRH of the spktrains
% 
load(mtffile, 'sprtmf', 'sprsmf', 'tmfaxis', 'smfaxis')

% CRH of the neurons
mtfhist = zeros(size(spktrain, 1), length(tmfaxis) * length(smfaxis));
if length(sprsmf) > length(spktrain)
    downsample = round(length(sprsmf)/length(spktrain));
    sprsmf = sprsmf(downsample:downsample:end);
    sprtmf = sprtmf(downsample:downsample:end);
end
diff = length(sprsmf) - length(spktrain);
if abs(diff) > 1
    error('stim subsample wrong!')
elseif diff > 0
    sprsmf = sprsmf(1:end-diff);
    sprtmf = sprtmf(1:end-diff);
elseif diff < 0
    spktrain = spktrain(:, 1:end+diff);
end

for ii = 1:size(spktrain, 1)
    tmf = rude(spktrain(ii,:), sprtmf);
    smf = rude(spktrain(ii,:), sprsmf);
    temp = histcounts2(smf, tmf, [smfaxis, 4], [tmfaxis, 64]);
    mtfhist(ii,:) = temp(:)';
end

