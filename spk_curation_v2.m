function spk_curation_v2(data_path, save_path)

% single unit curation for SqMon data from 2009-2013
% Congcong, 2022-01-21

%% get file names
listing = dir(data_path);
ListAll = struct2dataset(listing);

% delete the file names that do not have any data... (bytes<10000)
f = double(ListAll(:,4))>10000;
fileNames = ListAll(f,1);

for idxFile = 1:length(fileNames)
    %% get file info
    fileName = fileNames{idxFile,1};
    if ~contains(fileName, '.mat')
        continue
    end
    name_in = fullfile(data_path,fileName);
    load(name_in, 'spk', 'strf', 'trigger_dmr', 'trigger_dmrrep')
    fprintf('%s\n', name_in);
               
    %% make spk data structure from the pydict file
    
    spkoutfile = fullfile(save_path,sprintf('%s-curated.mat',fileName(1:end-4)));
    checkfile = dir(spkoutfile);
    if isempty(checkfile)
        for j = 1:length(spk)
            spk(j).spiketimes_dmr = double(spk(j).spiketimes_dmr) / spk(1).fs * 1000;
        end
        
        spk = calc_ISI_violation_percentage(spk, 1.5);
        spk = calc_firing_rate_stability(spk, trigger_dmr);
        
        % threshold spk based on ISI_vio and frs
        ISI_vio_thresh = 2;
        frs_thresh = 0.5;
        passed_unit = [spk.ISI_vio] < ISI_vio_thresh ...
            | [spk.frs] > frs_thresh;
        ndelete = length(passed_unit) - sum(passed_unit);
        if ndelete > 0
            fprintf('%d neurons do not pass the threshold\n', ndelete)
        end
            spk = spk(passed_unit);
        strf = strf(passed_unit);
        
        save(spkoutfile, 'spk', 'strf','trigger_dmr', 'trigger_dmrrep')
    end
end