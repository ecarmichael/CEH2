function shuffled_TC = MS_split_shuff(data_in, pos_in, keep_idx, nShuff, bins1, bins2)
%% MS_boot_shuff: run a bootstrap sampling for mean and CI of a data array.
%    Inputs:
%    - data_in: [1 x N array] data that will be bootstrapped to obtain mean
%    and CI
%
%    - pos_in:  [2 x N array]
%
%    - keep_idx [1 x N logical vector]  which indices to keep for analyses.
%     Typically this can be used to exclude moments of imobility.
%
%    - nShuff: number of shuffles.
%
%    - bins1: [vector]  vector of bins for dimension 1.
%
%    - bins2: [vector]  vector of bins for dimension 2 [optional]
%
%    Outputs:
%    - shuffled_TC [ n x nSamples]  shuffled tuning curve.  
%
%
% EC 2020-12-11   initial version  
% based off of Guillaume Etter's 2020 paper:
% https://www.frontiersin.org/articles/10.3389/fncir.2020.00019/full
%
%% initialize

if nargin == 1
    error('Requires data_in and at least one set of bins')
    
elseif nargin < 3
    keep_idx = ones(data_in); % keep all indicies if
    nShuff = 1000; % default number of shuffles.
    
elseif nargin < 4
    nShuff = 1000; % default number of shuffles.
elseif nargin < 5
    bins2 = []; % if only 1D then keep this empty.
end
%% use 1 or 2 methods.
% if isempty(bins2) % 1d
%     actual_bs_TC = zeros(length(bins1)-1, nShuff); % sub one to dims to account for centers.
%     shuffled_bs_TC = zeros(length(bins1)-1, nShuff);
% else
%     actual_bs_TC = zeros(length(bins1)-1,length(bins2)-1, nShuff);
%     shuffled_bs_TC = zeros(length(bins1)-1,length(bins2)-1, nShuff);
% end

% split shuffle the data.  (used for P values)  TODO: make this only for
% shuffled data.  and update for actual bootstrapping with random segments
% of the data.


for iShuff = nShuff:-1:1
    split_ts = ceil(MS_randn_range(1,1,1,length(data_in)));
    this_shuff = [data_in(end-split_ts+1:end); data_in(1:end-split_ts)]; % cut the data at a point and put the ends together.
    this_shuff_samples = ones(1,length(this_shuff));  % as a surrogate for timestamps.
    
    
    % get a set of indicies to include for this shuffle
    bootstrap_ts = ones(1,length(this_shuff));
    for ii = 1:length(bootstrap_ts)
        if bootstrap_ts(ii) == 1 && rand < 0.5
            bootstrap_ts(ii) = 0;
        end
    end
    bootstrap_ts = logical(bootstrap_ts);
    
    % Compute the actual tuning curve using a bootstrapped sample
    
    if isempty(bins2) % 1d methods
        
        [~,~,~,~, shuffled_TC(:,iShuff)]  = MS_get_spatial_information(data_in(bootstrap_ts),this_shuff_samples(bootstrap_ts), pos_in(bootstrap_ts,1), bins1);
        
    else % 2d methods.
        
        [~,~,~,~, shuffled_TC(:,:,iShuff)]  = MS_get_spatial_information_2D(this_shuff(bootstrap_ts),pos_in(bootstrap_ts,:), bins1,bins2);
        
    end
    
end % shuffled