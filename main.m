%% STEP0: spike curation---------------------------------------------------
datapath = '/data/congcong/SqMoPhys_Josh/mountainsort/pydict/std4_dmr'; %result of ms_output
savepath = '/data/congcong/SqMoPhys_Josh/mountainsort/pydict/std4_dmr'; 
trigpath = '/data/congcong/SqMoPhys_Josh/trigger';
cd(datapath)
% spk_curation(datapath, savepath, 'dmr', trigpath)
% save half ms binned spike trains
stimfolder = '/data/congcong/SqMoPhys_Josh/stim';
stimfile = fullfile(stimfolder, 'contra-dmr-176flo-20000fhi-4SM-64TM-15db-48DF-96khz-10min_DFt1_DFf5_stim.mat');
mtfstimfile = regexprep(stimfile, '_stim', '_mtf');
assert(~isempty(dir(mtfstimfile)), 'Get mtf file using batch_stimulus_to_tmf_smf first!')
load(mtfstimfile)
% files = dir('*-dmr*curated.mat');
% files = {files.name};
% badfiles = ne_batch_save_halfms_binned_spktrain(files, stimfile);

%% STEP1: plot results of sorting and response to dmr----------------------
cd('/data/congcong/SqMoPhys_Josh/mountainsort/pydict/std4_dmr')
%spkfiles = dir('*poly2*curated.mat');
spkfiles = dir('*-curated.mat');
figurepath = '/data/congcong/SqMoPhys_Josh/figure/singleunit/su_response_to_dmr_4std';
stimspr = fullfile(stimfolder, 'contra-dmr-176flo-20000fhi-4SM-64TM-15db-48DF-96khz-10min_DFt1_DFf5.spr');
plotflag = 0;
for ii = 1:length(spkfiles)
    clear waveform_STRF_CRH
    % load spike sorted result
    load(spkfiles(ii).name)
    if isempty(spk.spk)
        continue
    end
    if ~exist('waveform_STRF_CRH', 'var')
        % load the raw  data (filtered at 100-6000Hz )
        matfile = dir(sprintf('/data/congcong/SqMoPhys_Josh/multiunit/filt/*%s*', spk.exp));
        load(fullfile(matfile.folder, matfile.name), 'amplifier_data_filt');
        
        % get waveforms by aligning time stamps from sorted result and raw trace
        waveform = sm_get_waveforms(amplifier_data_filt, spk);
        
        % get strfs
        trigger = trigger(1,:);
        strf = calculate_strf(spk, trigger, 0, stimspr);
        mdb = strf(1).mdb;
        dur = (trigger(end)-trigger(1))/20000;
        taxis = strf(1).taxis;
        faxis = strf(1).faxis;
        
        % get ACG
        ACG = {};
        for jj = 1:length(spk.spk)
            acg = xcorr(spktrain(jj,:), 100); % at 0.5ms resolution, +-50ms
            acg(101) = 0;
            ACG{jj} = acg;
        end
        
        % get CRH
        [mtfhist, tmfedges, smfedges] = sm_calculate_CRH(spktrain, mtfstimfile);        
        
        % save waveform, strf and crh result in a structure
        waveform_STRF_CRH = strf;
        unit = {spk.spk.unit};
        [waveform_STRF_CRH.unit] =  unit{:};
        position = {spk.spk.position};
        [waveform_STRF_CRH.position] = position{:};
        avgwaveform = waveform(:,1);
        [waveform_STRF_CRH.avgwaveform] = avgwaveform{:};
        stdwaveform = waveform(:,2);
        [waveform_STRF_CRH.stdwaveform] = stdwaveform{:};
        [waveform_STRF_CRH.ACG] = ACG{:};
        for jj = 1:length(waveform_STRF_CRH)
            waveform_STRF_CRH(jj).mtfhist = mtfhist(jj,:);
            waveform_STRF_CRH(jj).tmfedges = tmfedges;
            waveform_STRF_CRH(jj).smfedges = smfedges;
            
        end
        save(spkfiles(ii).name, 'waveform_STRF_CRH', '-append')
    else
        
        strf = waveform_STRF_CRH;
        waveform = [{waveform_STRF_CRH.avgwaveform}; {waveform_STRF_CRH.stdwaveform}];
        waveform = waveform';
        ACG = {waveform_STRF_CRH.ACG};
        mtfhist = {waveform_STRF_CRH.mtfhist};
        mtfhist = cell2mat(mtfhist');
        taxis = strf(1).taxis;
        faxis = strf(1).faxis;
        tmfedges = strf(1).tmfedges;
        smfedges = strf(1).smfedges;
    end
    
    % plot properties of each neuron
    if contains(spk.probe, 'Tetrode') && plotflag
        % determine the location of waveforms in the plot
        x = 0:0.02:0.8;
        if strcmp(spk.probe, 'TetrodeB1x64')
            prbfile = '/data/congcong/SqMoPhys_Josh/prb/prbfiles/TetrodeB1x64.csv';
        elseif strcmp(spk.probe, 'Tetrode1x64')
            prbfile = '/data/congcong/SqMoPhys_Josh/prb/prbfiles/Tetrode1x64.csv';
        end
        prblayout = readmatrix(prbfile);
        plotx = prblayout(:,1)/100;
        idx = mod(plotx, 1) > 0.3;
        plotx(idx) = floor(plotx(idx))+ 1;
        idx = mod(plotx, 1) > 0.1;
        plotx(idx) = floor(plotx(idx))+ 0.5;
        ploty = (prblayout(:,2)-50)/75;
        idx = mod(ploty, 1) > 0.4;
        ploty(idx) = floor(ploty(idx))+ 1;
        idx = mod(ploty, 1) > 0.2;
        ploty(idx) = floor(ploty(idx))+ 0.5;
    elseif strcmp(spk.probe, 'poly21x48') && plotflag
        prbfile = '/data/congcong/SqMoPhys_Josh/prb/prbfiles/poly21x48.csv';
        prblayout = readmatrix(prbfile);
        prby =  prblayout(:,2);
        prby = sort(prby);
        ymid = prby(24);
        prblayout(prblayout(:,2) <= ymid, 1) = prblayout(prblayout(:,2) <= ymid, 1) + 2*86;
        prblayout(prblayout(:,2) > ymid, 2) = prblayout(prblayout(:,2) > ymid, 2) - ymid;
        plotx = prblayout(:,1)/86;
        ploty = prblayout(:,2)/50/4;
        prblayout = readmatrix(prbfile);
    end
    for jj = 1:length(spk.spk)
        figure('Renderer', 'painters', 'Position', [30 30 600 750]);
        % plot waveform of signal from each channel
        subplot(5,2,[1 3])
        avgwaveform = waveform{jj,1};
        avgwaveform = avgwaveform/max(max(abs(avgwaveform)))/2;
        hold on
        for kk = 1:length(prblayout)
            plot(x + plotx(kk), avgwaveform(kk,:) + ploty(kk))
        end
        position = abs(spk.spk(jj).position);
        positionidx = find(prblayout(:,1) == position(1) & prblayout(:,2) == position(2));
        positionx = plotx(positionidx);
        positiony = ploty(positionidx);
        text(positionx, positiony, '*', 'color', 'r')
        title(sprintf('%d spikes', spk.spk(jj).num_events))
        axis off
        hold off
        
        subplot(5,2,[2 4])
        bar(-50:0.5:50, ACG{jj}, 1, 'k')
        xticks([-50 50])
        title(sprintf('ISI violation: %.3f%%', spk.spk(jj).ISI_vio))
        ymax = max(ACG{jj}*1.1);
        if ymax < 1
            ymax = 1;
        end
        ylim([0 ymax])
        
        % plot STRF
        subplot(5,2,[5 7])
        sta = strf(jj).rfcontra;
        n0 = strf(jj).n0contra;
        plot_strf_raw(sta, faxis, taxis)
        if isfield(waveform_STRF_CRH, 'crh_sig_z')
            title(sprintf('STRF z:%d-p:%d', waveform_STRF_CRH(jj).strf_sig_z, waveform_STRF_CRH(jj).strf_sig_p(3)))
        else
            title('STRF')
        end
        
        subplot(5, 2, [6 8])
        crh = mtfhist(jj,:);
        plot_CRH(crh, tmfedges, smfedges)
        if isfield(waveform_STRF_CRH, 'crh_sig_z')
            title(sprintf('CRH z:%d-p:%d', waveform_STRF_CRH(jj).crh_sig_z, waveform_STRF_CRH(jj).crh_sig_p(3)))
        else
            title('CRH')
        end
        
        subplot(5,2, [9, 10])
        fr = spktrain(jj,1:2000*floor(length(spktrain)/2000));
        spknum1 = sum(fr(1:length(fr)/2));
        spknum2 = sum(fr(length(fr)/2+1:end));
        stability = (spknum2-spknum1)/(spknum1+spknum2);
        fr = reshape(fr, [2000, length(fr)/2000]);
        fr = sum(fr, 1);
        plot(fr)
        ylabel('Firirng rate (Hz)')
        xlabel('time (s)')
        title(sprintf('FR instability: %.3f', stability))
        waveform_STRF_CRH(jj).fr_stability = stability;
        saveas(gcf, fullfile(figurepath, regexprep(spkfiles(ii).name, '.mat', sprintf('-%d.jpg', jj))))
        close
        
    end
    
    save(spkfiles(ii).name, 'waveform_STRF_CRH', '-append')
end

%% STEP2: threshold units based on their firing rate stability------------
% plot_strf_significance_vs_firing_rate_stability;
spkfiles = dir('*-curated.mat');
for ii = 1:length(spkfiles)
    clear waveform_STRF_CRH
    load(spkfiles(ii).name, 'waveform_STRF_CRH')
    if exist('waveform_STRF_CRH', 'var')
        idx = find(abs([waveform_STRF_CRH.fr_stability]) < 0.5 & [waveform_STRF_CRH.w0contra] > 0.1);
        if isempty(idx)
            continue
        end
        load(spkfiles(ii).name, 'spk', 'trigger', 'spktrain', 'edges')
        waveform_STRF_CRH = waveform_STRF_CRH(idx);
        spk.spk = spk.spk(idx);
        spktrain = spktrain(idx,:);
        save(fullfile('/data/congcong/SqMoPhys_Josh/mountainsort/pydict/std4_dmr_thresh', [spkfiles(ii).name(1:end-4) '-thresh.mat']), 'trigger', 'spktrain', 'edges', 'spk', 'waveform_STRF_CRH')
        fprintf('%s\n', spkfiles(ii).name)
    end
end

%% STEP3: get the significance and properties for STRF -------------------
% calculate response properties of individual neurons that can be derived
% from STRF
% 1. get high resolution STRF using stim file that is not downsampled
% 2. calculate RI of strf
% 3. determine if there is robust feature of strf
cd('/data/congcong/SqMoPhys_Josh/mountainsort/pydict/std4_dmr_thresh')
spkfiles = dir('*thresh.mat');
stimfolder = '/data/congcong/SqMoPhys_Josh/stim';
sprfile = fullfile(stimfolder, 'contra-dmr-176flo-20000fhi-4SM-64TM-15db-48DF-96khz-10min.spr');
stimfile_downsampled = fullfile(stimfolder, 'contra-dmr-176flo-20000fhi-4SM-64TM-15db-48DF-96khz-10min_DFt1_DFf5_stim.mat');
if ~exist('stim', 'var')
    load(stimfile_downsampled, 'stim')
end

flag_plot = 1;
figurefolder_strf = '/data/congcong/SqMoPhys_Josh/figure/singleunit/su_response_to_dmr_4std_thresh/STRF';
figurefolder_RI = '/data/congcong/SqMoPhys_Josh/figure/singleunit/su_response_to_dmr_4std_thresh/RI';

nlags = 100;% for RI calculation, only 50ms preceding spikes are considered
niter = 1000;
parpool(2)
parfor ii = 1:length(spkfiles)
    data = load(spkfiles(ii).name, 'spk', 'waveform_STRF_CRH', 'trigger', 'spktrain')
%     if ~exist('waveform_STRF_CRH', 'var')
%         continue
%     end
    %spk = data.spk;
    waveform_STRF_CRH = data.waveform_STRF_CRH;
    spktrain = data.spktrain;
    fprintf('Processing %s\n', spkfiles(ii).name)
    % 1. calculate STRF without downsample
%     strf = calculate_strf(spk.spk, trigger(1,:), 0, sprfile);
%     tmp = {strf.spl}; [waveform_STRF_CRH.spl] = tmp{:};
%     tmp = {strf.taxis}; [waveform_STRF_CRH.taxis] = tmp{:};
%     tmp = {strf.faxis}; [waveform_STRF_CRH.faxis] = tmp{:};
%     tmp = {strf.rfcontra}; [waveform_STRF_CRH.rfcontra] = tmp{:};
%     tmp = {strf.spln}; [waveform_STRF_CRH.spln] = tmp{:};
    % 2. calculate reliability index of neurons
    waveform_STRF_CRH = calc_strf_RI_batch(waveform_STRF_CRH, spktrain, stim, nlags,...
        'plot_flag', 1, 'figure_folder', figurefolder_RI, 'figure_name_base', spkfiles(ii).name(1:13));
    
    % 3. for reliable STRFs, determine if the strf has robust features
%     fs = waveform_STRF_CRH.fs;
%     dur = (trigger(end) - trigger(1))/fs;
%     waveform_STRF_CRH = strf_properties(waveform_STRF_CRH, dur, flag_plot, ...
%         'save_plot', 1, 'figure_folder', figurefolder_strf, 'figure_name_base', spkfiles(ii).name(1:13));
    parsave(spkfiles(ii).name, waveform_STRF_CRH)
end

%% ---------------- step6: comparison of core and belt ------------------
rf_significance_core_vs_belt;
%% get distribution of firing rate instability
% spkfiles = dir('*-spk-curated.mat');
% fr_instability = [];
% for ii= 1:length(spkfiles)
%     clear waveform_STRF_CRH
%     load(spkfiles(ii).name', 'waveform_STRF_CRH')
%     if ~exist('waveform_STRF_CRH', 'var')
%         continue
%     end
%     
%     fr_instability = [fr_instability, [waveform_STRF_CRH.fr_stability]];
% end
% histogram(fr_instability, 50)
% xlabel('Firirng Rate Instability Index')
% ylabel('counts')
%% firing rate distribution
cd('/data/congcong/SqMoPhys_Josh/mountainsort/pydict/dmr_thresh')
spkfiles = dir('*poly2*.mat');
fr = [];
for ii = 1:length(spkfiles)
    clear waveform_STRF_CRH
    load(spkfiles(ii).name, 'waveform_STRF_CRH')
    if exist('waveform_STRF_CRH', 'var')
        if isfield(waveform_STRF_CRH, 'fr_stability')
            fr = [fr, [waveform_STRF_CRH.w0contra]];
        else
            fprintf('%s\n', spkfiles(ii).name)
        end
    end
end
%
histogram(fr, 0:0.4:20)
xlabel('Firing Rate')
ylabel('# of units')
title('poly2 1x48')

%% 
%spkfiles = dir('*-Tetrode*-curated.mat');
% numevents = [];
% isolation = [];
% noise_overlap = [];
% peak_snr = [];
% ISI_vio = [];
% for ii = 1:length(spkfiles)
%     load(spkfiles(ii).name)
%     numevents = [numevents [spk.spk.num_events]];
%     isolation = [isolation [spk.spk.isolation]];
%     noise_overlap = [noise_overlap [spk.spk.noise_overlap]];
%     peak_snr = [peak_snr [spk.spk.peak_snr]];
%     ISI_vio = [ISI_vio [spk.spk.ISI_vio]];
% end
%
% figure
% subplot(221)
% scatter(numevents, ISI_vio)
% title('ISI violation')
% subplot(222)
% scatter(numevents, isolation)
% title('isolation')
% subplot(223)
% scatter(numevents, peak_snr)
% title('peak SNR')
% subplot(224)
% scatter(numevents, noise_overlap)
% title('noise overlap')

%% ------ plot strf, rtf and crh--------------
cd('/data/congcong/SqMoPhys_Josh/mountainsort/pydict/std4_dmr_thresh')
spkfiles = dir('*-thresh.mat');
stimfolder = '/data/congcong/SqMoPhys_Josh/stim';
figurepath = '/data/congcong/SqMoPhys_Josh/figure/singleunit/su_response_to_dmr_4std_thresh';
stimspr = fullfile(stimfolder, 'contra-dmr-176flo-20000fhi-4SM-64TM-15db-48DF-96khz-10min_DFt1_DFf5.spr');
stimparamsfile = regexprep(stimspr, '.spr', '_param.mat');
load(stimparamsfile, 'faxis','MaxFM','MaxRD')
p = 0.02;
plotflag = 1;
for ii = 1:length(spkfiles)
    clear waveform_STRF_CRH
    % load spike sorted result
    load(spkfiles(ii).name)
    if isempty(spk.spk)
        continue
    end
    if ~exist('waveform_STRF_CRH', 'var')
        % load the raw  data (filtered at 100-6000Hz )
        matfile = dir(sprintf('/data/congcong/SqMoPhys_Josh/multiunit/filt/*%s*', spk.exp));
        load(fullfile(matfile.folder, matfile.name), 'amplifier_data_filt');
        
        % get waveforms by aligning time stamps from sorted result and raw trace
        waveform = sm_get_waveforms(amplifier_data_filt, spk);
        
        % get strfs
        trigger = trigger(1,:);
        strf = calculate_strf(spk, trigger, 0, stimspr);
        mdb = strf(1).mdb;
        dur = (trigger(end)-trigger(1))/20000;
        taxis = strf(1).taxis;
        faxis = strf(1).faxis;
        
        % get ACG
        ACG = {};
        for jj = 1:length(spk.spk)
            acg = xcorr(spktrain(jj,:), 100); % at 0.5ms resolution, +-50ms
            acg(101) = 0;
            ACG{jj} = acg;
        end
        
        % get CRH
        [mtfhist, tmfedges, smfedges] = sm_calculate_CRH(spktrain, mtfstimfile);        
        
        % save waveform, strf and crh result in a structure
        waveform_STRF_CRH = strf;
        unit = {spk.spk.unit};
        [waveform_STRF_CRH.unit] =  unit{:};
        position = {spk.spk.position};
        [waveform_STRF_CRH.position] = position{:};
        avgwaveform = waveform(:,1);
        [waveform_STRF_CRH.avgwaveform] = avgwaveform{:};
        stdwaveform = waveform(:,2);
        [waveform_STRF_CRH.stdwaveform] = stdwaveform{:};
        [waveform_STRF_CRH.ACG] = ACG{:};
        for jj = 1:length(waveform_STRF_CRH)
            waveform_STRF_CRH(jj).mtfhist = mtfhist(jj,:);
            waveform_STRF_CRH(jj).tmfedges = tmfedges;
            waveform_STRF_CRH(jj).smfedges = smfedges;
            
        end
        save(spkfiles(ii).name, 'waveform_STRF_CRH', '-append')
    else  
        strf = waveform_STRF_CRH;
        waveform = [{waveform_STRF_CRH.avgwaveform}; {waveform_STRF_CRH.stdwaveform}];
        waveform = waveform';
        ACG = {waveform_STRF_CRH.ACG};
        mtfhist = {waveform_STRF_CRH.mtfhist};
        mtfhist = cell2mat(mtfhist');
        taxis = strf(1).taxis;
        faxis = strf(1).faxis;
        tmfedges = strf(1).tmfedges;
        smfedges = strf(1).smfedges;
    end
    
    % plot properties of each neuron
    if contains(spk.probe, 'Tetrode') && plotflag
        % determine the location of waveforms in the plot
        x = 0:0.02:0.8;
        if strcmp(spk.probe, 'TetrodeB1x64')
            prbfile = '/data/congcong/SqMoPhys_Josh/prb/prbfiles/TetrodeB1x64.csv';
        elseif strcmp(spk.probe, 'Tetrode1x64')
            prbfile = '/data/congcong/SqMoPhys_Josh/prb/prbfiles/Tetrode1x64.csv';
        end
        prblayout = readmatrix(prbfile);
        plotx = prblayout(:,1)/100;
        idx = mod(plotx, 1) > 0.3;
        plotx(idx) = floor(plotx(idx))+ 1;
        idx = mod(plotx, 1) > 0.1;
        plotx(idx) = floor(plotx(idx))+ 0.5;
        ploty = (prblayout(:,2)-50)/75;
        idx = mod(ploty, 1) > 0.4;
        ploty(idx) = floor(ploty(idx))+ 1;
        idx = mod(ploty, 1) > 0.2;
        ploty(idx) = floor(ploty(idx))+ 0.5;
    elseif strcmp(spk.probe, 'poly21x48') && plotflag
        prbfile = '/data/congcong/SqMoPhys_Josh/prb/prbfiles/poly21x48.csv';
        prblayout = readmatrix(prbfile);
        prby =  prblayout(:,2);
        prby = sort(prby);
        ymid = prby(24);
        prblayout(prblayout(:,2) <= ymid, 1) = prblayout(prblayout(:,2) <= ymid, 1) + 2*86;
        prblayout(prblayout(:,2) > ymid, 2) = prblayout(prblayout(:,2) > ymid, 2) - ymid;
        plotx = prblayout(:,1)/86;
        ploty = prblayout(:,2)/50/4;
        prblayout = readmatrix(prbfile);
    end
    for jj = 1:length(spk.spk)
        figure
        figuresetup2savepdf(25, 25)
        % plot waveform of signal from each channel
        subplot(3,3,1)
        avgwaveform = waveform{jj,1};
        avgwaveform = avgwaveform/max(max(abs(avgwaveform)))/2;
        hold on
        for kk = 1:length(prblayout)
            plot(x + plotx(kk), avgwaveform(kk,:) + ploty(kk))
        end
        position = abs(spk.spk(jj).position);
        positionidx = find(prblayout(:,1) == position(1) & prblayout(:,2) == position(2));
        positionx = plotx(positionidx);
        positiony = ploty(positionidx);
        text(positionx, positiony, '*', 'color', 'r')
        title(sprintf('%d spikes', spk.spk(jj).num_events))
        axis off
        hold off
        
        subplot(3,3,2)
        bar(-50:0.5:50, ACG{jj}, 1, 'k')
        xticks([-50 50])
        title('ACG')
        ymax = max(ACG{jj}*1.1);
        if ymax < 1
            ymax = 1;
        end
        ylim([0 ymax])
        
        subplot(3,3,3)
        ISI = diff(spk.spk(jj).spiketimes);
        histogram(ISI, 0:1:100, 'FaceColor', 'k', 'FaceAlpha', 1)
        hold on
        plot([2 2], ylim, 'r--', 'linewidth', 2)
        title(sprintf('ISI violation: %.3f%%', spk.spk(jj).ISI_vio))

        % plot STRF
        subplot(3,3,4)
        sta = strf(jj).rfcontra;
        n0 = strf(jj).n0contra;
        plot_strf_raw(sta, faxis, taxis)
        if isfield(waveform_STRF_CRH, 'crh_sig_z')
            title(sprintf('STRF z:%d-p:%d', waveform_STRF_CRH(jj).strf_sig_z, waveform_STRF_CRH(jj).strf_sig_p(3)))
        else
            title('STRF')
        end
         
        % plot CRH
        subplot(3, 3, 6)
        crh = mtfhist(jj,:);
        tmfcenters = (tmfedges(2:end) + tmfedges(1:end-1))/2;
        smfcenters = (smfedges(2:end) + smfedges(1:end-1))/2;
        plot_CRH(crh, tmfcenters, smfcenters)
        if isfield(waveform_STRF_CRH, 'crh_sig_z')
            title(sprintf('CRH z:%d-p:%d', waveform_STRF_CRH(jj).crh_sig_z, waveform_STRF_CRH(jj).crh_sig_p(3)))
        else
            title('CRH')
        end
        
         % plot RTF
        subplot(3,3,5)
        nf = length(faxis);
        nlags = length(taxis);
        [tmf, xmf, rtf] = sm_mtf_sta2rtf(reshape(sta, nf, nlags), taxis, faxis, MaxFM, MaxRD);
        plot_CRH(fliplr(rtf), tmfedges, smfedges)
        title('rtf')
        
        subplot(3,3, 7:9)
        fr = spktrain(jj,1:2000*floor(length(spktrain)/2000));
        spknum1 = sum(fr(1:length(fr)/2));
        spknum2 = sum(fr(length(fr)/2+1:end));
        stability = (spknum2-spknum1)/(spknum1+spknum2);
        fr = reshape(fr, [2000, length(fr)/2000]);
        fr = sum(fr, 1);
        plot(fr)
        ylabel('Firirng rate (Hz)')
        xlabel('time (s)')
        title(sprintf('FR instability: %.3f', stability))
        waveform_STRF_CRH(jj).fr_stability = stability;
        saveas(gcf, fullfile(figurepath, regexprep(spkfiles(ii).name, '.mat', sprintf('-%d.jpg', jj))))
        close
        
    end
    
    save(spkfiles(ii).name, 'waveform_STRF_CRH', '-append')
end
