function [c, BF, BW_max, BW_min, Latency,Q, tw, sig, ex] = strf_bf_bw_latency(filter,taxis,faxis, TM, SM, varargin)
% calculate BF, BW, Latency for a filter (STA, MID1, MID)
% INPUTS
% filter: spectrotemporal receptive field (rfsig)
% taxis: time axis for the filter (s)
% faxis: frequency axis for the filter (Hz)
% bw_crit: the threshold for obtaining contour of response field
% OUTPUTS
% BF: best frequency (kHz)
%       list of values in order of:
%       [feature with the largest positive deviation;
%       feature with the largest negative deviation]
% BW_max: maximum frequency of the bandwidth (kHz)
% BW_min: minimum frequency of the bandwidth (kHz)
% Latency: latency of peak value (ms)
% BW: bandwidth (kHz), BW_max-BW_min
% Q: BF/BW
% tw: time window in which there is robust feature
% sig: if the strf has robust feature
% ex: 1 if the robust feature is excitatory, 2 if inhibitory

% Congcong, 2022-01-06

% if the SNR of the largest deviation is low, reject
if max(abs(filter(:))) < 3*std(filter(:))
    [c, BF, BW_max, BW_min, Latency, Q, tw, sig] = rejection;
    ex = NaN;
    return
end
% if the peak-trough difference before spike is not much higher than after,
% reject
% t0 = find(taxis == 0);
% if ~isempty(t0) && t0 > 1
%     nsample = t0-1;
%     pt_after = max(max(filter(:,1:nsample))) - min(min(filter(:,1:nsample)));
%     pt_before = max(max(filter(:,t0+1:t0+nsample))) - min(min(filter(:,t0+1:t0+nsample)));
%     SNR = pt_before/pt_after;
%     if SNR < 1.5
%         [c, BF, BW_max, BW_min, Latency, Q, tw, sig] = rejection;
%         ex = NaN;
%         return
%     end
% end

% if the maximum deviation is not within reasonable time bound (0-50ms), reject
[~, idx] = max(abs(filter(:)));
[~, col] = ind2sub(size(filter), idx);
delay = taxis(col);
if delay < 0 || delay > 0.05
    [c, BF, BW_max, BW_min, Latency, Q, tw, sig] = rejection;
    ex = NaN;
    return
end

p = inputParser;
% to make the width of the distribution at 1/e of the peak height (or 0.75, Atencio et al. 2012)
addParameter(p,'bw_crit',0.37);
parse(p,varargin{:});
bw_crit = p.Results.bw_crit;

% seperate the strf into excitatory and inhibitory subfields
pfilter = filter;
nfilter = filter;
pfilter(filter<0) = 0;
nfilter(filter>0) = 0;

% get properties of robust features in both excitatory and inhibitory
% fields
[ce, BFe, BW_maxe, BW_mine, Latencye,  Qe, twe, sige] = peak_contour(pfilter,faxis,taxis, TM, SM, bw_crit);
[cn, BFn, BW_maxn, BW_minn, Latencyn,  Qn, twn, sign] = peak_contour(nfilter,faxis,taxis, TM, SM, bw_crit);

if ~sige && ~sign
    [c, BF, BW_max, BW_min, Latency, Q, tw, sig] = rejection;
    ex = NaN;
else
    c = {ce, cn};
    BF = [BFe, BFn];
    BW_max = [BW_maxe; BW_maxn];
    BW_min = [BW_mine; BW_minn];
    Latency = [Latencye, Latencyn];
    Q = [Qe, Qn];
    tw = [twe; twn];
    sig = [sige, sign];
    [~, idx] = max(abs(filter));
    if filter(idx) > 0
        ex = 1;
    else
        ex = 2;
    end
end


end

function  [c, BF, BW_max, BW_min, Latency, Q, tw, sig] = rejection
    c = NaN;
    BF = NaN;
    BW_max = NaN;
    BW_min = NaN;
    Latency = NaN;
    Q = NaN;
    tw = NaN;
    sig = 0;
end

function  [c, BF, BW_max, BW_min, Latency, Q, tw, sig] = peak_contour(filter,faxis,taxis, TM, SM, bw_crit)
filter = abs(filter);
[peak, peak_idx] = max(filter(:));
[row, col] = ind2sub(size(filter), peak_idx);
Latency = taxis(col) * 1000;
if Latency < 0 || Latency > 50
    [c, BF, BW_max, BW_min, Latency, Q, tw, sig] = rejection;
    return
end
filter_binary = filter;
filter_binary(filter > 0) = 1;
M = contourc(filter_binary, 1);
idx = 1;
while idx < size(M, 2)
    nvertices = M(2, idx);
    c = M(:, idx+1:idx+nvertices);
    [in, on] = inpolygon(col, row, c(1,:), c(2,:));
    if in || on
        break
    end
    idx = idx + nvertices + 1; 
end
% how large is the feature on time axis
difft = taxis(2) - taxis(1);
tw = (max(c(1,:))-min(c(1,:))) *  difft * 1000; %in ms
t_thresh = 1/TM/2*0.76 * 1000;
if tw < t_thresh
    [c, BF, BW_max, BW_min, Latency, Q, tw, sig] = rejection;
    return
end
% how large is the feature on frequency axis
difff = log2(faxis(end) / faxis(1))/ (length(faxis)-1);
fw = (max(c(2,:))-min(c(2,:))) *  difff; %in octave
f_thresh = 1/SM/2*0.76;
if fw < f_thresh
    [c, BF, BW_max, BW_min, Latency, Q, tw, sig] = rejection;
    return
end

% if the size of the feature pass the threshold
filter(filter >= peak*bw_crit) = 1;
filter(filter < peak*bw_crit) = 0;
M = contourc(filter, 1);
idx = 1;
while idx < size(M, 2)
    nvertices = M(2, idx);
    c = M(:, idx+1:idx+nvertices);
    [in, on] = inpolygon(col, row, c(1,:), c(2,:));
    if in || on
        break
    end
    idx = idx + nvertices + 1; 
end

BF = faxis(row)/1000;
BW_max = faxis(round(max(c(2,:))))/1000;
BW_min = faxis(round(min(c(2,:))))/1000;
Q = (BW_max - BW_min)/BF;
tw = (max(c(1,:)) - min(c(1,:)))*1000;
sig = 1;
end