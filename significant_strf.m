function [rfsig] = significant_strf(rf, p, n0, mdb, dur, soundtype)
% [rfsig] = significant_strf(rf, p, n0, mdb, dur, soundtype)
%
% Get significant strf at the p value for dynamic moving ripple
% stimulus.
%
% rf is the non-significant strf.
%
% p is the significance level.
%
% n0 is the number of spikes used to make the strf.
%
% mdb is the modulation depth of the ripple stimulus.
%
% dur is the duration of the ripple stimulus in sec, 
%    and this may be computed as 
%       (trigger(end)-trigger(1) ) / fs, 
%    where fs is the sampling rate of the A/D system.
%    The default value is 15 min * 60 sec = 900 sec
%
% soundtype is either 'dmr*' for dynamic moving ripple or
%    'rn*' for ripple noise. The '*' means wildcard.
%
% If dur is not given then 15 minutes will be assumed.
%
% Formerly named wstrfstat.m and located in Monty's Keck Toolbox.
%
% caa 7/7/03


if ( nargin < 4 | nargin > 6 )
   error('significant_strf.m needs 4 to 6 inputs args.');
end

if ( nargin == 4 )
   dur = 15 * 60;
   soundtype = 'dmr';
end

if ( nargin == 5 )
   soundtype = 'dmr';
end


if ( ~isempty(findstr(soundtype, 'dmr')) )
   sigma2 = (mdb)^2 / 8; % variance of dynamic moving ripple
else
   sigma2 = (mdb)^2 / 12; % variance of dynamic moving ripple
end


sigma = sqrt(sigma2); % std of dynamic moving ripple

siglevel = norminv(1-p/2) * sqrt(n0) / ( sigma * dur );

rfabs = abs(rf);

rfsig = zeros(size(rf));

rfsig(rfabs > siglevel) = rf(rfabs > siglevel);



