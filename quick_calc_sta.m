function [sta, stabigmat, spkcountvec] = quick_calc_sta(stimulus, locator, varargin)

% Quick calculate STA (using matrix multiplication) based on the system's
% memory.

% Inputs:
%   stimulus: (nf * nt) stimulus matrix.
%   locator: (m * nt) spike train matrix, where m = number of neurons.
%   nlags: number of time bins to look at before each spike.
%   chunks: number of chunks to calculate STA with

% Outputs:
%   sta: sta for each neuron (1 per row)
%   stabigmat: sta for each chunk of stimulus/spktrain segment (neurons x
%   sta size x number of chunks)
%   spkcountvec: spike count for each chunk of spike train

% Written by JS, 12/7/16.
% Updated by JS, 6/18/18, included 'chunks' option to manually input
% number of chunks. Included stabigmat as an output.
% Updated by JS, 3/18/19, to suppress output when required.
% Updated by JS, 10/07/19, to include 'nleads'
% Updated by CH, 01/14/2022, to use GPU computing

ip = inputParser;
addRequired(ip, 'stimulus', @ismatrix);
addRequired(ip, 'locator', @(x) ismatrix(x) || isvector(x));
addParameter(ip, 'nlags', 20, @(x) mod(x,1) == 0);
addParameter(ip, 'chunks', [], @(x) mod(x,1) == 0);
addParameter(ip, 'suppressprint', 0, @(x) x == 0 || x == 1);
addParameter(ip, 'nleads', 2, @(x) mod(x,1) == 0);
addParameter(ip, 'GPU', 0, @(x) x == 0 || x == 1);

parse(ip, stimulus, locator, varargin{:});

fn = fieldnames(ip.Results);

for i = 1:length(fn)
    eval(sprintf('%s = ip.Results.%s;', fn{i}, fn{i}));
end

assert(size(stimulus,2) == size(locator,2));

nf = size(stimulus,1);

if isempty(chunks)
    % calculate max array based on memory
    if GPU
        GPUparameter = gpuDevice();
        user.MaxPossibleArrayBytes = GPUparameter.AvailableMemory;
    else
        try
            user = memory;
        catch
            user.MaxPossibleArrayBytes = 3e+10;
        end
    end
    MaxElements = user.MaxPossibleArrayBytes / 8;
    limit = MaxElements / 5;
    
    % dimensions of stim matrix
    size1 = size(stimulus,2) - nlags + 1;
    size2 = nlags * nf;
    numelements = size1*size2;
    
    % determine number of chunks
    chunks = ceil(numelements / limit);
    
end

chunksize = floor(size(stimulus,2) / chunks);

% get chunk indices
startidx = [1 chunksize.*(1:chunks-1)-nlags + 2]; % include previous nlag bins for calculation of STA
endidx = [chunksize.*(1:chunks-1) size(stimulus,2)];
if GPU
    locator = gpuArray(locator);
    
    stabigmat = zeros(size(locator,1), nf * (nlags+nleads), chunks, 'gpuArray');
    spkcountvec = zeros(size(locator, 1), chunks, 'gpuArray');
    
    for i = 1:chunks
        [stim, resp] = create_stim_trial_from_stim_matrix(stimulus(:,startidx(i):endidx(i)),...
            locator(:,startidx(i):endidx(i)), nlags, nleads); 
        spkcountvec(:,i) = sum(locator(:,startidx(i):endidx(i)), 2);
        stabigmat(:,:,i) = resp*stim;
    end
    
else
    stabigmat = zeros(size(locator,1), nf * (nlags+nleads), chunks);
    spkcountvec = zeros(size(locator, 1), chunks);
    
    if ~suppressprint
        fprintf('\nNumber of chunks for STA calculation: %d', chunks)
    end
    
    for i = 1:chunks
        
        if ~suppressprint
            fprintf('\nCalculating chunk %d of %d', i, chunks)
        end
        [stim, resp] = create_stim_trial_from_stim_matrix(stimulus(:,startidx(i):endidx(i)),...
            locator(:,startidx(i):endidx(i)), nlags, nleads);
        
        spkcountvec(:,i) = sum(locator(:,startidx(i):endidx(i)), 2);
        stabigmat(:,:,i) = resp*stim;
        
    end
end
sta = sum(stabigmat, 3);

if ~suppressprint
    fprintf('\n')
end

