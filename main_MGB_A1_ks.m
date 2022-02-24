%% STEP0: spike curation---------------------------------------------------
spk_path = '/data/congcong/rat_MGB_A1/3_singleunit/dmr';
cd(spk_path)
% save half ms binned spike trains
files = dir('*-newspk.mat');
files = {files.name};
%stimfolder = '/data/congcong/SqMonPhys_Brian/stim';
%stimfile = fullfile(stimfolder, 'contra-dmr-176flo-20000fhi-4SM-64TM-15db-48DF-96khz-10min_DFt1_DFf5_stim.mat');
%badfiles = ne_batch_save_halfms_binned_spktrain(files);
%ne_batch_save_halfms_binned_spktrain_for_spon_evoked_activity(files)
%% STEP1: get waveform of each neuron----------------------
cd('/data/congcong/rat_MGB_A1/3_singleunit/dmr')
spkfiles = dir('*-newspk.mat');
figfolder_waveform = '/data/congcong/rat_MGB_A1/figure/singleunit/waveform';
ks2_output_path = '/data/congcong/rat_MGB_A1/2_spikesorted/kilosort/output';

flag_plot = 1;
for ii = 1:length(spkfiles)
    
    v = whos('-file', spkfiles(ii).name);
    if 0%sum(strcmp({v.name}, 'waveform')) == 1
        fprintf('(%d/%d) Waveform already exist %s\n', ii, length(spkfiles), spkfiles(ii).name)
        continue
    end
    % load spike sorted result
    data = load(spkfiles(ii).name, 'newspk', 'spktrain');
    spk = data.newspk;
    spktrain = data.spktrain;
    
    if isempty(spk.spk)
        fprintf('(%d/%d) No neuron in %s\n', ii, length(spkfiles), spkfiles(ii).name)
        continue
    end
    fprintf('(%d/%d) Processing %s\n', ii, length(spkfiles), spkfiles(ii).name)
    % load the raw  data (filtered at 100-6000Hz )
    binfile = dir(sprintf('/data/congcong/rat_MGB_A1/2_spikesorted/kilosort/input/%s*%s*/%s*-filtered.bin', spk.exp, spk.probe, spk.exp));
    bytes = get_file_size(fullfile(binfile.folder, binfile.name));
    nTimepoints = floor(bytes/64/2); % number of total timepoints
    fid = fopen(fullfile(binfile.folder, binfile.name));
    fullArray = fread(fid,[64 nTimepoints],'int16');
    fclose(fid);
    
    % get ks channel infromation
    KSdir = fullfile(ks2_output_path, spkfiles(ii).name(1:end-35));
    params.cgsFile = fullfile(KSdir, 'cluster_KSLabel.tsv');
    params.excludeNoise = false;
    ss = loadKSdir(KSdir, params);
    if length(ss.chs) < 64
        %keyboard
        fullArray = fullArray(ss.chs,:);
    end
    % get waveforms by aligning time stamps from sorted result and raw trace
    waveform_avg_std = sm_get_waveforms(fullArray, spk);
    
    % get ACG(during dmr)
    ACG = cell(length(spk.spk), 1);
    for jj = 1:length(spk.spk)
        acg = xcorr(spktrain(jj,:), 100); % at 0.5ms resolution, +-50ms
        acg(101) = 0;
        ACG{jj} = acg;
    end
    
    n0contra = num2cell(cellfun(@(x) length(x), {spk.spk.spiketimes_dmr}));
    % get basic information from spk
    waveform = struct('exp', {spk.exp}, ...
                      'unit', {spk.spk.unit},...
                      'chan', {spk.spk.chan}, ...
                      'position', {spk.spk.position}, ...
                      'n0contra', n0contra, ...
                      'w0contra', {spk.spk.fr_dmr});
    avgwaveform = waveform_avg_std(:,1);
    [waveform.avgwaveform] = avgwaveform{:};
    stdwaveform = waveform_avg_std(:,2);
    [waveform.stdwaveform] = stdwaveform{:};
    taxis = waveform_avg_std(:,3);
    [waveform.taxis] = taxis{:};
    [waveform.ACG] = ACG{:};
    
    basename = fullfile(figfolder_waveform, spkfiles(ii).name(1:13));
    waveform = waveform_properties(waveform, 'flag_plot', flag_plot,...
        'fig_basename', basename, 'index0', 0);
    % get firing rate stability and ISI
    for jj = 1:size(spktrain,1)
        % fr stability
        halfsample = floor(size(spktrain,2)/2);
        spknum1 = sum(spktrain(jj, 1:halfsample));
        spknum2 = sum(spktrain(jj, halfsample+1:end));
        stability = (spknum2-spknum1)/(spknum1+spknum2);
        
        % ISI
        ISI = diff(spk.spk(jj).spiketimes);
        ISI_vio = sum(ISI < 1.5)/waveform(jj).n0contra;
        
        waveform(jj).fr_stability = stability;
        waveform(jj).ISI = ISI;
        waveform(jj).ISI_vio = ISI_vio;
    end    

    parsave(spkfiles(ii).name, 'waveform', waveform)
end

%% STEP3: get the significance and properties for STRF -------------------
% calculate response properties of individual neurons that can be derived
% from STRF
% 1. get high resolution STRF using stim file that is not downsampled
% 2. calculate RI of strf
% 3. determine if there is robust feature of strf
cd('/data/congcong/rat_MGB_A1/3_singleunit/dmr')
spkfiles = dir('*newspk.mat');
stimfolder = '/data/congcong/stimulus/thalamus';
sprfile = fullfile(stimfolder, 'rn1-500flo-40000fhi-0-4SM-0-40TM-40db-96khz-48DF-15min-seed190506.spr');
stimfile_downsampled = fullfile(stimfolder, 'rn1-500flo-40000fhi-0-4SM-0-40TM-40db-96khz-48DF-15min-seed190506_DFt1_DFf5_stim.mat');
if ~exist('stim_mat', 'var')
    load(stimfile_downsampled, 'stim_mat')
end
figurefolder_strf = '/data/congcong/rat_MGB_A1/figure/singleunit/STRF';
figurefolder_RI = '/data/congcong/rat_MGB_A1/figure/singleunit/STRF/RI';
nlags = 100;% for RI calculation, only 50ms preceding spikes are considered
parfor ii = 1:length(spkfiles)
    data = load(spkfiles(ii).name, 'newspk', 'waveform', 'trigger', 'spktrain');
    waveform = data.waveform;
    spktrain = data.spktrain;
    trigger = data.trigger;
    spk = data.newspk;
    %strf = data.strf;
    fprintf('Processing %s\n', spkfiles(ii).name)
    % 1. calculate STRF without downsample
    strf = calculate_strf(spk.spk, trigger(1,:), 0, sprfile);
    % 2. calculate reliability index of neurons
%     strf = calc_strf_RI_batch(strf, spktrain, stim, nlags,...
%         'plot_flag', 1, 'figure_folder', figurefolder_RI, 'figure_name_base', spkfiles(ii).name(1:13));
%     
    % 3. for reliable STRFs, determine if the strf has robust features
    fs = spk.fs;
    dur = (trigger(end) - trigger(1))/fs;
    strf = strf_properties(strf, dur, 'plot_flag', 1, ...
        'save_plot', 1, 'figure_folder', figurefolder_strf, 'figure_name_base', spkfiles(ii).name(1:13));
    parsave(spkfiles(ii).name, 'strf', strf)
end

%% STEP4: get the properties for rtf -------------------
cd('/data/congcong/rat_MGB_A1/3_singleunit/dmr')
spkfiles = dir('*newspk.mat');

flag_plot = 1;
figurefolder_rtf = '/data/congcong/rat_MGB_A1/figure/singleunit/rtf';

for ii = 1:length(spkfiles)
    % load data
    fprintf('Processing %s\n', spkfiles(ii).name)
    data = load(spkfiles(ii).name, 'strf', 'trigger');
    strf = data.strf;
    trigger = data.trigger;
    % get rtf
    fs = strf(1).fs;
    dur = (trigger(end) - trigger(1))/fs;
    basename = fullfile(figurefolder_rtf, sprintf('%s', spkfiles(ii).name(1:13)));
    rtfparams = strf_parameters(strf, dur, flag_plot, basename);
    
    parsave(spkfiles(ii).name, 'rtfparams', rtfparams)
end

%% STEP5: get the significance and properties for CRH -------------------
% get mtf file
% cd('/data/congcong/SqMoPhys_Josh/stim')
% batch_stimulus_to_tmf_smf({'contra-dmr-176flo-20000fhi-4SM-64TM-15db-48DF-96khz-10min_DFt1_DFf5_stim.mat'},
% 100, 8, 16); % 50ms seconds window before spike
cd('/data/congcong/rat_MGB_A1/3_singleunit/dmr')
spkfiles = dir('*newspk.mat');
stimfolder = '/data/congcong/stimulus/thalamus';
mtffile = fullfile(stimfolder, 'rn1-500flo-40000fhi-0-4SM-0-40TM-40db-96khz-48DF-15min-seed190506_mtf.mat');
load(mtffile, 'sprtmf', 'sprsmf', 'tmfaxis', 'smfaxis')

flag_plot = 1;
figurefolder_crh = '/data/congcong/rat_MGB_A1/figure/singleunit/crh';
figurefolder_RI = '/data/congcong/rat_MGB_A1/figure/singleunit/crh/RI';
parfor ii = 1:length(spkfiles)
    
    fprintf('Processing %s\n', spkfiles(ii).name)
    data = load(spkfiles(ii).name, 'newspk', 'spktrain');
    spk = data.newspk;
    spktrain = data.spktrain;
    
    % get crh
    [mtfhist, tmfaxis, smfaxis] = sm_calculate_CRH(spktrain, mtffile);
    
    % save parameters to crh structure
    mtfhist = mat2cell(mtfhist, ones(1, size(mtfhist, 1)), size(mtfhist, 2));
    crh = struct('exp', {spk.exp}, ...
                 'unit', {spk.spk.unit},...
                 'chan', {spk.spk.chan}, ...
                 'position', {spk.spk.position}, ...
                 'tmfaxis', tmfaxis, ...
                 'smfaxis', smfaxis, ...
                 'mtfhist', mtfhist');
    % get crh properties
    basename = fullfile(figurefolder_crh, sprintf('%s', spkfiles(ii).name(1:13)));
    crh = batch_crh_parameters(crh, flag_plot, basename);
    % get crh RI
    basename = fullfile(figurefolder_RI, sprintf('%s', spkfiles(ii).name(1:13)));
%     crh = batch_crh_RI(crh, sprtmf, sprsmf, spktrain,...
%         'plot_flag', 1, 'fig_basename', basename, 'nblocks', 10);
    
    parsave(spkfiles(ii).name, 'crh', crh)
end

%% ------ plot strf, rtf and crh--------------
cd('/data/congcong/SqMoPhys_Josh/mountainsort/pydict/std4_dmr_thresh')
spkfiles = dir('*-thresh.mat');
figurefolder = '/data/congcong/SqMoPhys_Josh/figure/singleunit';
for ii = 1:length(spkfiles)
    % load spike sorted result
    data = load(spkfiles(ii).name,'spk', 'waveform', 'strf', 'rtfparams', 'crh', 'trigger');
    waveform = data.waveform;
    strf = data.strf;
    rtfparams = data.rtfparams;
    crh = data.crh;
    probe = data.spk.probe;
    trigger = data.trigger;
    basename = fullfile(figurefolder, sprintf('%s', spkfiles(ii).name(1:13)));
    plot_su_strf_rtf_crh(probe, waveform, strf, rtfparams, crh, trigger, basename)
end
