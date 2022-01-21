function pli = phase_locking_index(strf, mdb, w0, riptype)
% phase_locking_index - calculate phase locking index for strfs.
%
% strf : significant portion of strf found using the
%    function significant_strf.m, though not the rate-
%    normalized strf.
%
% mdb : modulation depth of ripple stimulus used to
%    obtain the strf.
%
% n0 : number of spikes used to construct the strf.
%
% w0 : mean spike rate, spikes per second.
%
% riptype : ripple stimulus type. A string that may
%    be 'dmr*' or 'rn*' where * is wildcard.
%
% caa 7/7/03


if ( nargin ~= 4 )
   error('You need 4 input arguments.');
end


if ( ~isempty(findstr(riptype, 'dmr')) )
   sigma = mdb / sqrt(8);
   delta = sqrt(8);
elseif ( ~isempty(findstr(riptype,'rn')) )
   sigma = mdb / sqrt(12);
   delta = sqrt(12);
else
   error('riptype must be dmr or rn.');
end


pli = sigma*(max(max(strf)) - min(min(strf))) / (w0 * delta);

return;

