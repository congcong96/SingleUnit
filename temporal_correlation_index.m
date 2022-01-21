function [r, tshift] = temporal_correlation_index(rf, delta)
% temporal_correlation_index - temporal correlation in strf.
%
% r = temporal_correlation_index(rf, fs)
%
% r is between 0 and 1, where 0 implies the strf 
% in rf has a completely random temporal structure
% and 1 implies the rf is temporally uniform. r is
% calculated for shifts of 5, 10, 15, ..., 50 ms.
%
% tshift : number of msec bin shifts for which r is calculated.
%    Will be 1x10 vector, where each element of shift
%    represents a 5, 10, 15, ..., or 50 ms shift.
%
% rf : strf matrix
%
% delta : difference in time between time samples of strf
%
% The temporal_correlation_index is from 
% Touryan, Lau, Dan (2002). Isolation of relevant 
% features from random stimuli for cortical
% complex cells. J Neurosci 22(24):10811-8
%
% caa 7/10/03


tshift = 0:.001:0.05; % 5 ms bin shifts; largest shift is 100 ms
r = zeros(size(tshift));
for i = 1:length(tshift)
   bins = ceil(tshift(i) / delta);
   rfb = circshift(rf,bins ,2);  % shift temporally, not spectrally
   r(i) = corr2(rf,rfb);  % calculate the correlation value
end

return;


