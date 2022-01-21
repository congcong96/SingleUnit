function [energy] = strf_energy(rf, mdb, riptype)
%
% [energy] = strf_energy(rf)
%
% FILE NAME       : STRF CORR COEF
% DESCRIPTION     : Computes the energy in the strf. The
%    statistically significant, though not the rate-normalized
%    strf, is assumed to be the input.
%	
% rf		: Spectro-temporal receptive field
%
% energy      : Standard deviation of STRF1
%

if ( ~isempty(findstr(riptype, 'dmr')) )
   sigma = mdb / sqrt(8);
   delta = sqrt(8);
elseif ( ~isempty(findstr(riptype,'rn')) )
   sigma = mdb / sqrt(12);
   delta = sqrt(12);
else
   error('riptype must be dmr or rn.');
end

rf = sigma * rf; % rate-normalized strf

%Finding Non-Zero values and rearanging STRFs to a linear array
index = find( rf ~= 0 );
X1 = rf(index);

% Finding Mean Energy.
energy = sqrt(mean(X1.^2));

return;

