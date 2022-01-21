function rel_idx = calc_CRH_RI_v0(sprtmf, sprsmf, spktrain)

% calculate reliability index based on JS paper2
% crh is first calculated on 10 chunks of data, and then 5 of the chunks
% are added to correlate with the sum of the other 5 chunks
% the 5 chunks are selected randomly and the radom sampling is done 100
% time and the 100 correlation values are later averaged to get RI
% Inputs:
%   sprtmf: temporal modulation frequency of stimulus
%   sprsmf: spectral modulation frequency of stimulus
%   spktrain: locators of spike train (or cNE event train) with the same
%   length and temporal resolution as sprtmf
% output:
%   rel_idx: 100 correlation values for each unit
% Congcong, 2021-05-24

numchunk = 10; % number of chunks on which to calculate crh
numbins = floor(size(spktrain,2)/numchunk);
for ii = 1:numchunk
    spktrain_tmp = spktrain(:, (ii-1)*numbins+1 : ii*numbins);
    sprtmf_tmp = sprtmf((ii-1)*numbins+1 : ii*numbins);
    sprsmf_tmp = sprsmf( (ii-1)*numbins+1 : ii*numbins);
    crh(ii,:, :) = calc_crh(spktrain_tmp, sprtmf_tmp, sprsmf_tmp);
end

rel_idx = zeros(size(spktrain,1), 100);
for ii = 1:100
    
    idxA = randsample(1:numchunk, numchunk/2);
    idxB = setdiff(1:numchunk, idxA);
    
    assert(length(idxA) == length(idxB))
    
    crhA = squeeze(sum(crh(idxA, :, :), 1));
    crhB = squeeze(sum(crh(idxB, :, :), 1));
    [ro, col] = size(crhA);
    if ro < col
        crhA = crhA'; 
        crhB = crhB';
    end
    
    rel_idx(:,ii) = diag(corr(crhA, crhB));
end

end

function crh = calc_crh(spktrain, sprtmf, sprsmf)

tmf_hist = -64:4:64;
smf_hist = 0:.2:4;

for ii = 1:size(spktrain,1)
    tmf = rude(spktrain(ii,:), sprtmf);
    smf = rude(spktrain(ii,:), sprsmf);
    temp = histcounts2(smf, tmf, smf_hist, tmf_hist);
    crh(ii,:) = temp(:);
end

end