function [taxis,faxis,STRF1,STRF2,PP,Wo1,Wo2,No1,No2,SPLN,SIG1,SIG2]=rtwstrfdb(SpecFile,T1,T2,spet,Trig,Fss,SPL,MdB,ModType,Sound,NBlocks,sprtype)

%function [taxis,faxis,STRF1,STRF2,PP,Wo1,Wo2,No1,No2,SPLN]=rtwstrfdb(SpecFile,
%T1,T2,spet,Trig,Fss,SPL,MdB,ModType,Sound,NBlocks,sprtype)
%
%       FILE NAME       : RT WSTRF DB
%       DESCRIPTION     : Real Time spectro-temporal receptive field
%			  Uses Lee/Schetzen Aproach via Specto-Temporal Envelope
%			  For dB Amplitude Sound distributions 
%
%	SpecFile	: Spectral Profile File
%	T1, T2		: Evaluation delay interval for WSTRF(T,F)
%			  T E [- T1 , T2 ], Note that T1 and T2 > 0
%	spet		: Array of spike event times in sample number
%	Trig		: Array of Trigger Times
%	Fss		: Sampling Rate for TRIGGER and SPET
%	SPL		: Signal RMS Sound Pressure Level
%	MdB		: Signal Modulation Index in dB
%	ModType		: Kernel modulation type : 'lin' or 'dB'
%	Sound		: Sound Type 
%			  Moving Ripple	: MR ( Default )
%			  Ripple Noise	: RN
%	NBlocks		: Number of Blocks Between Displays
%	sprtype		: SPR File Type : 'float' or 'int16'
%			  Default=='float'	
%
%	RETURNED VALUES 
%
%	taxis		: Time Axis
%	faxis		: Frequency Axis (Hz)
%	STRF1 , STRF2	: Spectro-Temporal Receptive Field
%	PP		: Power Level
%	Wo1, Wo2	: Zeroth-Order Kernels ( Average Number of Spikes / Sec )
%	No1, No2	: Number of Spikes
%	SPLN		: Sound Pressure Level per Frequency Band
%
% add a calculation of robustness of strf based on the difference of responses between minimum and maximum values for stimullus window compared to
% spontaneous window (see details below).
% Natsumi 30Aug17

%Parameters
if nargin<12
    sprtype='float';
end
NTrig=length(Trig);

%Loading Parameter Data
index=strfind(SpecFile,'.spr');
ParamFile=[SpecFile(1:index(1)-1) '_param.mat'];
load(ParamFile, 'Fs', 'DF', 'NF', 'NT', 'faxis')

%Fliping Trig and Spet for channel 2 and channel 1 STRFs
spet = spet(:)';
Trig = Trig(:)';
MinTime=min([Trig spet]);
MaxTime=max([Trig spet]);
spet=spet-MinTime+1;
Trig=Trig-MinTime+1;

%Converting Temporal Delays to Sample Numbers
N1=round(T1*Fs/DF); % DF = temporal downsampling factor for spectral profile 
N2=round(T2*Fs/DF);
 
%Opening Spectral Profile File
fid=fopen(SpecFile);

%Initializing Some Variables
TrigCount=2;
No1=0;				%Number of Spikes for channel 1
No2=0;				%Number of Spikes for channel 2
STRF1=zeros(NF,N1+N2);		%Receptive Field  for channel 1
STRF2=zeros(NF,N1+N2);		%Receptive Field  for channel 2

%Computing Spectro Temporal Receptive Fields - Checking for 'dB' or 'lin'
if strcmp(ModType,'dB')

	%Fiding Mean Spectral Profile and RMS Power
	SPLN=SPL-10*log10(NF);					% Normalized SPL per frequency band
	if strcmp(Sound,'RN')
		RMSP=-MdB/2;					% RMS value of normalized Spectral Profile % if MdB = 40, the modulation will be in the range from -20 dB to 20 dB.
		PP=MdB^2/12;					% Modulation Depth Variance 
	elseif strcmp(Sound,'MR')
		RMSP=-MdB/2;					% RMS value of normalized Spectral Profile
		PP=MdB^2/8;					% Modulation Depth Variance 
	end

	%Initializing First and Second Spectral Profile Segments
	frewind(fid);
	if strcmp(sprtype,'float')
		S1=fread(fid,NT*NF,'float');
		S2=fread(fid,NT*NF,'float');
		S3=fread(fid,NT*NF,'float');
	else
		S1=fread(fid,NT*NF,'int16')/.99/1024/32/2-.5;
		S2=fread(fid,NT*NF,'int16')/.99/1024/32/2-.5;
		S3=fread(fid,NT*NF,'int16')/.99/1024/32/2-.5;
	end
	S1=reshape(S1,NF,NT);
	S2=reshape(S2,NF,NT);
	S3=reshape(S2,NF,NT);

	%Loading Data and Computing 'dB' STRF
	while ~feof(fid) && TrigCount<length(Trig)-1

		%Finding SPET in between triggers
		index1=spet>=Trig(TrigCount) & spet<Trig(TrigCount+1);
		index2=spet>Trig(NTrig-TrigCount+1) & spet<=Trig(NTrig-TrigCount+2);

		%Resampling spet relative to the Spectral Profile samples
		spettrig1=ceil( (spet(index1)-Trig(TrigCount)+1) * Fs / Fss /DF );
		spettrig2=ceil( (Trig(NTrig-TrigCount+2)+1-spet(index2)) * Fs / Fss /DF );

		%Finding Receptive Field for Channel 1
		epsilon=10^(-MdB/20);
		for k=1:length(spettrig1)

			%Setting Spike Time and STRF length
			M=size(S1);M=M(2);
			L=spettrig1(k);

			%Averaging Pre-Event Spectral Profiles
			if L < N2
				STRF1=STRF1+ MdB*[S1(:,M-(N2-L-1):M) S2(:,1:L+N1)] - RMSP;
			elseif L+N1 > M
				STRF1=STRF1+ MdB*[S2(:,L-N2+1:M) S3(:,1:N1-M+L)] - RMSP;
			else
				STRF1=STRF1+ MdB*S2(:,L-N2+1:L+N1) - RMSP;
			end

			%Counting the number of Spikes averaged
			No1=No1+1;

		end

		%Finding Receptive Field for Channel 2
		for k=1:length(spettrig2)
	
			%Setting Spike Time and STRF length
			M=size(S1);M=M(2);
			L=spettrig2(k);

			%Averaging Pre-Event Spectral Profiles
			if L < N1
				STRF2=STRF2+ MdB*[S1(:,M-(N1-L-1):M) S2(:,1:L+N2)] - RMSP;
			elseif L+N2 > M
                try
  				STRF2=STRF2+ MdB*[S2(:,L-N1+1:M) S3(:,1:N2-M+L)] -RMSP;
                catch
                    warning('Bad spet2');
                end
			else
				STRF2=STRF2+ MdB*S2(:,L-N1+1:L+N2) - RMSP;
			end

			%Counting the number of Spikes averaged
			No2=No2+1;

		end

		%Reading Spectral Profile Data File
		S1=S2;
		S2=S3;
		if strcmp(sprtype,'float')
			S3=fread(fid,NT*NF,'float');
		else
			S3=fread(fid,NT*NF,'int16')/.99/1024/32/2-.5;
		end
		if ~feof(fid)
			S3=reshape(S3,NF,NT);
		end

		%Updating Trigger Counter
		TrigCount=TrigCount+1;
% 		clc
      if ( ~mod(TrigCount, 100) )
   		disp(['Block Number ' num2str(TrigCount) ' of ' num2str(NTrig)])
      end 

% 		%Sending To Display - Updates every NBlocks
		if TrigCount/NBlocks==round(TrigCount/NBlocks)
			T=min([ ( Trig(TrigCount) - Trig(1) )/Fss  (max(spet) - min(spet))/Fss]);
			Wo1=No1/T;
			Wo2=No2/T;  
			taxis=(-N1:N2-1)/(Fs/DF);
% 
% 
% 			%Displaying STRF1
%          clf;
% 			subplot(2,1,1);
% 			pcolor(taxis,log2(faxis/500),Wo1/PP*fliplr(STRF1)/No1*sqrt(PP));
% 			shading flat;
%          cmap = brewmaps('rdbu', 21);
%          colormap(cmap);
% %          colormap jet,colorbar
% 			title(['No = ' int2str(No1) ' ( Spikes ), Wo = ' num2str(Wo1,5) ' ( Spikes/Sec )'])
% 
% 			%Displaying STRF2
% 			subplot(2,1,2);
% 			pcolor(taxis,log2(faxis/500),Wo2/PP*STRF2/No2*sqrt(PP));
% 			shading flat;
%          cmap = brewmaps('rdbu', 21);
%          colormap(cmap);
% % 			shading flat,colormap jet,colorbar
% 			title(['No = ' int2str(No2) ' ( Spikes ) , Wo = ' num2str(Wo2,5) ' ( Spikes/Sec )'])
% 			pause(0)
		end
	end

	%Normalizing 'dB' STRF According to Paper by Van Dijk
	if max(spet)>max(Trig)
		T=( max(Trig) - min(Trig) )/Fss;
    else  
		T=( max(spet) - min(spet) )/Fss;
	end
	if ~( T==0 || isempty(T) )
		Wo1=No1/T;
		Wo2=No2/T;
	else
		No1=-9999;						%No Spikes in SPET
		No2=-9999;
	end
	STRF1=Wo1/PP*fliplr(STRF1)/No1;%summation of STRF/PP/T?
	STRF2=Wo2/PP*STRF2/No2;
    taxis=(-N1:N2-1)/(Fs/DF);

elseif strcmp(ModType,'lin')

	%Fiding Mean Spectral Profile and RMS Power
        Po=2.2E-5;							% Threshold of Hearing at 1KHz in Pascals
        P= Po*10^(SPL/20);						% Pressure conversion
        PP=P*P/NF;							% Mean Power spectrum per frequency band
	SPLN=SPL-10*log10(NF);						% Normalized SPL per frequency band
	if strcmp(Sound,'RN')
%		a=10^(-MdB/10);
%		RMSP=(a-1)/log(a);                                      % RMS value of normalized Spectral Profile
%		MeanP=mean(10.^(MdB*fread(fid,1024*256,'float')/10));   % Same as RMSP but need to compute Numerically 
									% because Theoretical Slightly off
		a=10^(-MdB/20);
		RMSP=sqrt((a^2-1)/2/log(a)-((a-1)/log(a))^2);		% RMS value of normalized Spectral Profile
		MeanP=mean(10.^(MdB*fread(fid,1024*1024,'float')/20));  % Same as RMSP but need to compute Numerically 
									% because Theoretical Slightly off
	elseif strcmp(Sound,'MR')
		%dx=.01;
		%RMSP=1/2/pi*10^(-MdB/20);
		%RMSP=RMSP*dx*sum(10.^(MdB/20*sin(0:dx:2*pi)));          % RMS value of normalized Spectral Profile
		%MeanP=RMSP;
        
		RMSP=std(10.^(MdB*fread(fid,1024*1024,'float')/20));
		MeanP=mean(10.^(MdB*fread(fid,1024*1024,'float')/20));
	end
        
	%Initializing First and Second Spectral Profile Segments
	frewind(fid);
	if strcmp(sprtype,'float')
		S1=fread(fid,NT*NF,'float');
		S2=fread(fid,NT*NF,'float');
		S3=fread(fid,NT*NF,'float');
	else
		S1=fread(fid,NT*NF,'int16')/.99/1024/32/2-.5;
		S2=fread(fid,NT*NF,'int16')/.99/1024/32/2-.5;
		S3=fread(fid,NT*NF,'int16')/.99/1024/32/2-.5;
	end
	S1=reshape(S1,NF,NT);
	S2=reshape(S2,NF,NT);
	S3=reshape(S2,NF,NT);

	%Loading Data and Computing 'lin' STRF
	while ~feof(fid) && TrigCount<length(Trig)

		%Finding SPET in between triggers
		index1=spet>=Trig(TrigCount) & spet<Trig(TrigCount+1);
		index2=spet>Trig(NTrig-TrigCount+1) & spet<=Trig(NTrig-TrigCount+2);

		%Resampling spet relative to the Spectral Profile samples
		spettrig1=ceil( (spet(index1)-Trig(TrigCount)+1) * Fs / Fss /DF );
		spettrig2=ceil( (Trig(NTrig-TrigCount+2)+1-spet(index2)) * Fs / Fss /DF );

		%Finding Receptive Field for Channel 1
		for k=1:length(spettrig1)
	
			%Setting Spike Time and STRF length
			M=size(S1);M=M(2);
			L=spettrig1(k);

			%Averaging Pre-Event Spectral Profiles
			if L < N2
				%STRF1=STRF1+ ( 10.^(MdB*[S1(:,M-(N2-L-1):M) S2(:,1:L+N1)]/10) - MeanP ) / RMSP*PP;
				STRF1=STRF1+ ( 10.^(MdB*[S1(:,M-(N2-L-1):M) S2(:,1:L+N1)]/20) - MeanP ) / RMSP*sqrt(PP);
			elseif L+N1 > M
				%STRF1=STRF1+ ( 10.^(MdB*[S2(:,L-N2+1:M) S3(:,1:N1-M+L)]/10) - MeanP ) / RMSP*PP;
				STRF1=STRF1+ ( 10.^(MdB*[S2(:,L-N2+1:M) S3(:,1:N1-M+L)]/20) - MeanP ) / RMSP*sqrt(PP);
			else
				%STRF1=STRF1+ ( 10.^(MdB*S2(:,L-N2+1:L+N1)/10) - MeanP ) / RMSP*PP;
				STRF1=STRF1+ ( 10.^(MdB*S2(:,L-N2+1:L+N1)/20) - MeanP ) / RMSP*sqrt(PP);
			end

			%Counting the number of Spikes averaged
			No1=No1+1;

		end

		%Finding Receptive Field for Channel 2
		for k=1:length(spettrig2)

			%Setting Spike Time and STRF length
			M=size(S1);M=M(2);
			L=spettrig2(k);

			%Averaging Pre-Event Spectral Profiles
			if L < N1
				%STRF2=STRF2+ ( 10.^(MdB*[S1(:,M-(N1-L-1):M) S2(:,1:L+N2)]/10) - MeanP ) / RMSP*PP;
				STRF2=STRF2+ ( 10.^(MdB*[S1(:,M-(N1-L-1):M) S2(:,1:L+N2)]/20) - MeanP ) / RMSP*sqrt(PP);
			elseif L+N2 > M
				%STRF2=STRF2+ ( 10.^(MdB*[S2(:,L-N1+1:M) S3(:,1:N2-M+L)]/10) - MeanP ) / RMSP*PP;
				STRF2=STRF2+ ( 10.^(MdB*[S2(:,L-N1+1:M) S3(:,1:N2-M+L)]/20) - MeanP ) / RMSP*sqrt(PP);
			else
				%STRF2=STRF2+ ( 10.^(MdB*S2(:,L-N1+1:L+N2)/10) - MeanP ) / RMSP*PP;
				STRF2=STRF2+ ( 10.^(MdB*S2(:,L-N1+1:L+N2)/20) - MeanP ) / RMSP*sqrt(PP);
			end

			%Counting the number of Spikes averaged
			No2=No2+1;

		end

		%Reading Spectral Profile Data File
		S1=S2;
		S2=S3;
		if strcmp(sprtype,'float')
			S3=fread(fid,NT*NF,'float');
		else
			S3=fread(fid,NT*NF,'int16')/.99/1024/32/2-.5;
		end
		if ~feof(fid)
			S3=reshape(S3,NF,NT);
		end

		%Updating Trigger Counter
		TrigCount=TrigCount+1;
		clc
		disp(['Block Number ' num2str(TrigCount) ' of ' num2str(NTrig)])

		%Sending To Display - Updates every NBlocks
		if TrigCount/NBlocks==round(TrigCount/NBlocks)
	
			T=min([ ( Trig(TrigCount) - Trig(1) )/Fss  (max(spet) - min(spet))/Fss]);
			Wo1=No1/T;
			Wo2=No2/T;
			taxis=(-N1:N2-1)/(Fs/DF);
% 
% 			%Displaying STRF1
% 			subplot(211)
% 			%pcolor(taxis,log2(faxis/500),Wo1/fact(2)/PP^2*fliplr(STRF1)/No1*PP)
% 			pcolor(taxis,log2(faxis/500),Wo1/PP*fliplr(STRF1)/No1*sqrt(PP))
% 			shading flat,colormap jet,colorbar
% 			title(['No = ' int2str(No1) ' ( Spikes ) , Wo = ' num2str(Wo1,5) ' ( Spikes/Sec )'])
% 
% 			%Displaying STRF2
% 			subplot(212)
% 			%pcolor(taxis,log2(faxis/500),Wo2/fact(2)/PP^2*STRF2/No2*PP)
% 			pcolor(taxis,log2(faxis/500),Wo2/PP*STRF2/No2*sqrt(PP))
% 			shading flat,colormap jet,colorbar
% 			title(['No = ' int2str(No2) ' ( Spikes ) , Wo = ' num2str(Wo2,5) ' ( Spikes/Sec )'])
% 			pause(0)
		end
	end


	%Normalizing STRF According to Paper by Van Dijk
        if max(spet)>max(Trig)
                T=( max(Trig) - min(Trig) )/Fss;
        else
                T=( max(spet) - min(spet) )/Fss;
        end
	if ~( T==0 || isempty(T) )
		Wo1=No1/T;
		Wo2=No2/T;
	else
		No1=-9999;						%No Spikes
		No2=-9999;
	end
	%STRF1=Wo1/fact(2)/PP^2*fliplr(STRF1)/No1;
	%STRF2=Wo2/fact(2)/PP^2*STRF2/No2;
	STRF1=Wo1/PP*fliplr(STRF1)/No1;
	STRF2=Wo2/PP*STRF2/No2;
    taxis=(-N1:N2-1)/(Fs/DF);

end

% Calculation of robustness of strf based on the criteria: 1.5 times larger
% response for stimullus window compared to spontaneous window.

% find max and min values for spontanous and stimulus windows
tSpon = find(taxis<0);
tStim = find(taxis>=0);

for idx_window = 1:2
    switch idx_window
        case 1
            ftime = tSpon;
        case 2
            ftime = tStim;
    end
    Diffs1(idx_window) = max(max(STRF1(:,ftime))) -  min(min(STRF1(:,ftime)));
    Diffs2(idx_window) = max(max(STRF2(:,ftime))) -  min(min(STRF2(:,ftime)));
end
SIG1 = Diffs1(2)/Diffs1(1);
SIG2 = Diffs2(2)/Diffs2(1);

%Closing all opened files
fclose all;
