function newspk = ms_curate_data_based_on_params(spk, varargin)

% suggested values:
% isolation: 0.9
% noise_overlap: 0.05
% peak_snr: 2
% ISI_vio: 2
% numspikes: 100
% reject_thresh: 1


ip = inputParser;
addRequired(ip, 'spk', @isstruct)
addParameter(ip, 'isolation', [], @(x) isempty(x) || (x >= 0 && x <= 1))
addParameter(ip, 'noise_overlap', [], @(x) isempty(x) || (x >= 0 && x <= 1))
addParameter(ip, 'peak_snr', [], @(x) isempty(x) || isscalar(x))
addParameter(ip, 'firing_rate', [], @(x) isempty(x) || isscalar(x))
addParameter(ip, 'ISI_vio', [], @(x) isempty(x) || isscalar(x))
addParameter(ip, 'peak_amp', [], @(x) isempty(x) || isscalar(x))
addParameter(ip, 'reject_thresh', [], @(x) isempty(x) || isscalar(x))
addParameter(ip, 'numspikes', [], @(x) isempty(x) || isscalar(x))
parse(ip, spk, varargin{:})

spk = ip.Results.spk;
isolation = ip.Results.isolation;
noise_overlap = ip.Results.noise_overlap;
peak_snr = ip.Results.peak_snr;
firing_rate = ip.Results.firing_rate;
ISI_vio = ip.Results.ISI_vio;
peak_amp = ip.Results.peak_amp;
numspikes = ip.Results.numspikes;
reject_thresh = ip.Results.reject_thresh;

paramcount = 0;

allthresh = [];

if ~isempty(isolation)
    isovals = [spk.spk.isolation];
    isovalthresh = isovals < isolation;
    allthresh = [allthresh; isovalthresh];
    paramcount = paramcount + 1;
    filt_params.isolation = isolation;
end

if ~isempty(peak_snr)
    psnrvals = [spk.spk.peak_snr];
    psnrvalthresh = psnrvals < peak_snr;
    allthresh = [allthresh; psnrvalthresh];
    paramcount = paramcount + 1;
    filt_params.peak_snr = peak_snr;
end

if ~isempty(noise_overlap)
    nolvals = [spk.spk.noise_overlap];
    nolvalthresh = nolvals > noise_overlap;
    allthresh = [allthresh; nolvalthresh];
    paramcount = paramcount + 1;
    filt_params.noise_overlap = noise_overlap;
end

if ~isempty(firing_rate)
    frvals = [spk.spk.firing_rate];
    frvalthresh = frvals > firing_rate;
    allthresh = [allthresh; frvalthresh];
    paramcount = paramcount + 1;
    filt_params.firing_rate = firing_rate;
end

if ~isempty(ISI_vio)
    isivvals = [spk.spk.ISI_vio];
    isivvalthresh = isivvals > ISI_vio;
    allthresh = [allthresh; isivvalthresh];
    paramcount = paramcount + 1;
    filt_params.ISI_vio = ISI_vio;
end

if ~isempty(peak_amp)
    pavals = [spk.spk.peak_amp];
    pavalthresh = pavals < peak_amp;
    allthresh = [allthresh; pavalthresh];
    paramcount = paramcount + 1;
    filt_params.peak_amp = peak_amp;
end

if ~isempty(numspikes)
    nsvals = [spk.spk.num_events];
    frvalthresh = nsvals < numspikes;
    allthresh = [allthresh; frvalthresh];
    paramcount = paramcount + 1;
    filt_params.firing_rate = numspikes;
end
    

if isempty(reject_thresh)
    reject_thresh = floor(0.5*paramcount);
    filt_params.reject_thresh = reject_thresh;
end

spk.filt_params = filt_params;
    
rejectedunits = sum(allthresh,1) >= reject_thresh;
newspk = spk;
newspk.spk(rejectedunits) = [];


