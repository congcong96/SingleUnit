function strf = calculate_strf(spk, trigger, binaural, envfile, tbefore, tafter)
% 
% Simple script to call the strf function strfdbcalc.m
%
% strf = calculate_strf(spk, trigger, binaural, envfile, tbefore, tafter); 
%
% spk : struct array holding spike time data for the strf
%    calculations.
%
% trigger : array holding trigger times for the strf 
%    calculations.
%
% binaural : optional. Specifies if the stimulus was delivered
%            binaurally. 1 = yes, 0 = no. Default is 0.
%
% envfile : optional argument. Specifies the envelope file
%             to use for reverse correlation. The default is
%             #4, listed below.
%
% tbefore and tafter are optional arguments. They specify how much
% time you wish to process before a spike and how much time after.
% They must be in seconds. Example tafter = 0.050 and tbefore = 0.200,
% which means calculate an strf using stimulus data 50 ms before a
% spike occurs and 200 ms before a spike occurred. The default
% is tbefore = 0.200 and tafter = 0.050.
%
% caa 6/28/02

% --updated 11/07/2019 by Congcong
% ifchar(atten) -> str2double
% strf -> 1x64 struct

if ( nargin < 5 )
   T1 = 0.050; % 50 ms 
   T2 = 0.200; % 200 ms
elseif ( nargin == 6 )
   T1 = tafter;
   T2 = tbefore;
else
   error('You need either 4 or 6 input args.');
end

i2 = strfind(envfile,'SM');
i3 = strfind(envfile,'TM');
i4 = strfind(envfile,'db');

sm = str2double(envfile(i2(end)-1));
tm = str2double(envfile(i3(end)-2:i3(end)-1));
mdb = str2double(envfile(i4(end)-2:i4(end)-1));

if isfield(spk, 'atten')
    atten = spk.atten;
else
    atten = 0;
end
if ischar(atten)
    atten = str2double(atten(1:end-2));
end
if isfield(spk, 'fs')
    fs = double(spk(1).fs);
else
    fs = 20000;
end
spl = 105 - atten;

try
    stype = spk(1).stim;
    if iscell(spk(1).stim)
        stype = 'dmr';
    end
catch
    stype = 'dmr';
end

if ( contains(stype,'dmr') )
   snd = 'MR';
elseif ( contains(stype,'rn') )
   snd = 'RN';
else
   error('stim in struct array "spk" must be "dmr" or "rn".');
end

modtype = 'dB';
nblocks = 500;  % update display every 500 blocks

if ~isfield(spk, 'spiketimes')
%     strf = spk;
%     strf = rmfield(strf, 'spk');
    spk = spk.spk;
end


for i = 1:length(spk)
    fprintf('Processing neuron #%d\n', i)
   sp = spk(i).spiketimes;
   spet = round(sp / 1000 * double(fs)); % convert ms to sample number

   [taxis, faxis, rf1, rf2, pp, w01, w02, n01, n02, spln] = ...
      rtwstrfdb(envfile, T1, T2, spet', trigger, fs, spl, mdb, modtype, snd, nblocks);
  %STRF = summation of STRF/pp/T, see Monty 2001, JNeuro, in units of spk/s


   if ( binaural )

%       strf(i).exp = spk(i).exp;
%       strf(i).site = spk(i).site;
      strf(i).chan = spk(i).chan;
%       strf(i).model = spk(i).model;
%       strf(i).stim = spk(i).stim;
%       strf(i).atten = spk(i).atten;
%       strf(i).depth = spk(i).depth;
      if isfield(spk, 'position')
          strf(i).position = spk(i).position;
      end
      strf(i).sm = sm;
      strf(i).tm = tm;
      strf(i).mdb = mdb;
      strf(i).spl = spl;
      strf(i).taxis = taxis;
      strf(i).faxis = faxis;
      strf(i).n0contra = n01;
      strf(i).w0contra = w01;
      strf(i).rfcontra = single(rf1); % use single precision to save memory
      strf(i).n0ipsi = n02;
      strf(i).w0ipsi = w02;
      strf(i).rfipsi = single(rf2);
      strf(i).pp = pp;
      strf(i).spln = spln;
%       strf.strf(i).fs = fs;

   else  % it's monaural data
      strf(i).chan = spk(i).chan;

      if isfield(spk, 'position')
          strf(i).position = spk(i).position;
      end
      strf(i).sm = sm;
      strf(i).tm = tm;
      strf(i).mdb = mdb;
      strf(i).spl = spl;
      strf(i).taxis = taxis;
      strf(i).faxis = faxis;
      strf(i).n0contra = n01;
      strf(i).w0contra = w01;
      strf(i).rfcontra = single(rf1); % use single precision to save memory
      strf(i).pp = pp;
      strf(i).spln = spln;
      strf(i).fs = fs;

   end

   %fprintf('%d of %d channels completed.\n', i, length(spk));
   %pause(0.5);

end % (for)


