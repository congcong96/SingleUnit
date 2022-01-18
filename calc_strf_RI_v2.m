function [sig_z, sig_p, RI, z, p, RI_thresh, strfcorr, nstrfcorr] = calc_strf_RI_v2(stimulus, spktrain, nlags, varargin)
% Determine if the STRF is stable over time
% The RI of the STRF is compared to the null distribution of RIs of nSTRFs
% STRFs are calculated using quick_calc_sta
% null correlation distribution is generated by flipping the spiketrain
%
% Inputs:
%   stimulus: nf * nt matrix of stimulus
%   spktrain: spike train of neurons. (binned as stim, one at a time)
%   nlags: number of datapoints on taxis for strf calculation
%   niter: number of iteration to calculate reliability index, 
%             default = 1000
%   
% Outputs:
%   sig_z: significance of the reliability index of the unit based on
%   z-score (according to p 0.01, two-tail)
%   sig_p: significance based on p-value (0.95, 0.99, 0.999)
%   RI: reliability index of the unit
%   strfcorr: (RI_iter * 1) distribution of correlation of the neuron
%   nstrfcorr: (RI_iter * 1) null distribution of correlation
%
% Written by Congcong, 05/06/2021.
    

% RELIABILITY INDEX OF STRF
p = inputParser;
addParameter(p,'niter', 1000);
addParameter(p,'plot_flag', 0);
addParameter(p,'figure_folder', 0);
addParameter(p,'figure_name', 0);
parse(p,varargin{:});
niter = p.Results.niter;
plot_flag = p.Results.plot_flag;
figure_folder = p.Results.figure_folder;
figure_name = p.Results.figure_name;

spktrain1 = zeros(niter, length(spktrain));
spktrain2 = zeros(niter, length(spktrain));
idx = find(spktrain > 0);
n = floor(length(idx)/2);
for ii = 1:niter
    
    randidx = idx(randperm(length(idx)));
    
    spktrain1_tmp = spktrain;
    spktrain1_tmp(randidx(1:n))=0;
    spktrain1(ii,:) = spktrain1_tmp;
    
    spktrain2_tmp = spktrain;
    spktrain2_tmp(randidx(n+1:end))=0;
    spktrain2(ii,:) = spktrain2_tmp;
end
sta1 = quick_calc_sta(stimulus, spktrain1, 'chunks', 1, 'nlags', nlags, 'suppressprint', 1);
sta2 = quick_calc_sta(stimulus, spktrain2, 'chunks', 1, 'nlags', nlags, 'suppressprint', 1);
strfcorr = diag(corr(sta1', sta2'));
RI = mean(strfcorr);

% null distribution of flipped spiketrain
sta1 = quick_calc_sta(fliplr(stimulus), spktrain1, 'chunks', 1, 'nlags', nlags, 'suppressprint', 1);
sta2 = quick_calc_sta(fliplr(stimulus), spktrain2, 'chunks', 1, 'nlags', nlags, 'suppressprint', 1);
nstrfcorr = diag(corr(sta1', sta2'));

RI_thresh = [quantile(nstrfcorr, 0.95), quantile(nstrfcorr, 0.99), quantile(nstrfcorr, 0.999)];
p = sum(nstrfcorr >= RI)/niter;
sig_p =  RI > RI_thresh;

z = (mean(strfcorr) - mean(nstrfcorr))/sqrt((std(strfcorr)^2 + std(strfcorr)^2)/2);
sig_z =  z > 2.33;

if plot_flag
    figure
    histogram(strfcorr, 25, 'normalization', 'probability')
    hold on
    histogram(nstrfcorr, 25, 'normalization', 'probability')
    legend({'original', 'flipped'})
    xlabel('RI')
    ylabel('proportion')
    saveas(gcf, fullfile(figure_folder, figure_name))
    close
end





