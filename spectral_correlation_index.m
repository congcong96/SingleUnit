function [rx, xshift] = spectral_correlation_index(rf, deltax)
% spectral_correlation_index - spectral correlation in strf.
%
% [rx, xshift] = spectral_correlation_index(rf, deltax)
%
% rx is between 0 and 1, where 0 implies the strf 
% in rf has a completely random spectral structure
% and 1 implies the rf is spectrally uniform. r is
% calculated for different octave shifts.
%
% xshift : number of octave bin shifts for which r is calculated.
%    Will be 1x16 vector, where each element of shift
%    represents a -2.0, -1.75, ..., 1.75, or 2.0 octave shift.
%    A shift of zero is not included.
%
% rf : strf matrix
%
% deltax : difference, in octaves, between frequency samples of strf
%
% The spectral_correlation_index is adapted from 
% Touryan, Lau, Dan (2002). Isolation of relevant 
% features from random stimuli for cortical
% complex cells. J Neurosci 22(24):10811-8
%
% caa 7/10/03

% Positive shift means a frequency row is shifted
% to a higher frequency.

% xshift = [-2.0 -1.75 -1.5 -1.25 -1.0 -0.75 -0.5 -0.25 0 ...
% 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0];
xshift = -2:0.01:2;

for i = 1:length(xshift)
   bins = ceil(xshift(i) / deltax);
   rfb = circshift(rf,bins);  % shift spectrally, not temporally
   rx(i) = corr2(rf,rfb);  % calculate the correlation value
end

xshift = -xshift; % give output so that positive shifts
                  % correspond to upward shifts in frequency


return;





