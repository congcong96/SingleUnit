function [tmf, xmf, rtf] = sm_mtf_sta2rtf(v, taxis, faxis, MaxFM, MaxRD, tmodbins, smodbins, Display)
% sm_mtf_sta2rtf  Ripple transfer function from an STA, or other filter
% 
%     [tmf, xmf, rtf] = sm_mtf_sta2rtf(v, taxis, faxis, MaxFM, MaxRD, Display)
% 
%     Calculates the ripple transfer function from a filter.
% 
%     v : filter. May be either STA, MID1, or MID2. 
%     taxis   : Time Axis
%     faxis   : Frequency Axis
%     MaxFm   : Maximum Modulation Rate for Experiment
%     MaxRD   : Maximum Ripple Density for Experiment
%     Display : Display : 'y' or 'n'
%               Default : 'n'
% 
%     tmf : temporal modulation frequency axis (cycles / s).
%     xmf : spectral modulaiton frequency axis (cycles / octave).
%     rtf : ripple transfer function. A matrix where spectral modulation
%        varies along the rows and temporal modulation varies along the
%        columns of the matrix.



if ( nargin < 6 )
    tmodbins = 16;
    smodbins = 10;
   Display = 'n';
elseif  ( nargin < 8 )
   Display = 'n';
end

% Sampling rate along temporal/spectral filter axes
FsT = 1/(taxis(3)-taxis(2));
FsX = 1/log2(faxis(2)/faxis(1));

dtf = MaxFM / tmodbins; % Want 16 bins over the modulation range
ntbins = ceil(FsT / dtf); % sampling rate over temp res
ntbins = ntbins + ~rem(ntbins,2); % add a bin if we have an even number
%ntbins = 129;

dff = MaxRD / smodbins; % Want 20 bins over the modulation range
nfbins = ceil(FsX / dff); % sampling rate over freq resolution
nfbins = nfbins + ~rem(nfbins,2); % add a bin to get an odd number
%nfbins = 129;



% Estimate rtf using 2D FFT
rtf = fft2(v, nfbins, ntbins);

% Shift so 0 frequency is in the middle
rtf = fftshift( real( rtf.*conj(rtf) ) );


% Get the tm frequency vector - will be in cycles per second
if ( mod(size(rtf,2),2) )
   tmf = ( -(size(rtf,2)-1)/2:(size(rtf,2)-1)/2 ) / size(rtf,2) * FsT;
else
   tmf = ( -size(rtf,2)/2:(size(rtf,2)/2-1) ) / size(rtf,2) * FsT;
end


% Get the sm frequency vector - will be in cycles per octave
if ( mod(size(rtf,1),2) )
   xmf = [-(size(rtf,1)-1)/2:(size(rtf,1)-1)/2]/size(rtf,1) * FsX;
else
   xmf = [-size(rtf,1)/2:(size(rtf,1)/2-1)]/size(rtf,1) * FsX;
end

%Discarding Unecessary Data Samples
index1 = find( xmf >= 0 & xmf <= MaxRD );
index2 = find(tmf <= MaxFM & tmf >= -MaxFM);
rtf = rtf(index1,index2);
xmf = xmf(index1);
tmf = tmf(index2);


if ( strcmp(Display,'y') )
   figure; 
   subplot(1,2,1)
   imagesc(taxis, faxis, v)
   axis xy;
   subplot(1,2,2)
   imagesc(tmf, xmf, rtf )
   axis xy;
   set(gca,'tickdir', 'out');
   pause
   close all
end


return;

