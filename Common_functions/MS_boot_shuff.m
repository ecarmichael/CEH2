function stats_out = MS_boot_shuff(data_in, pos_in, keep_idx, nShuff, bins1, bins2)
%% MS_boot_shuff: run a bootstrap sampling to get the CI of a data array (can be 1D or 2D).  shuffles the data and takes 50% of the samples (repeated nShuff times)
%
%
%
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
%    - stats_out: [struct]  contains the cfg struct as well as stats values
%    for mean, SD, SEM, CI.
%
%
%
%
% EC 2020-12-11   initial version  
% based off of Guillaume Etter's 2020  paper: https://www.frontiersin.org/articles/10.3389/fncir.2020.00019/full
%
%
%
%% initialize

if nargin == 1
    error('Requires data_in and at least one set of bins')
    
elseif nargin < 3
    keep_idx = ones(data_in); % keep all indicies if
    nShuff = 1000; % default number of shuffles.
    
elseif nargin < 4
    nShuff = 1000; % default number of shuffles.
elseif nargin < 6
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
    
    
    % Get 50% of the indicies for this shuffle.
    
    bootstrap_ts = zeros(length(data_in),1);
    bootstrap_ts(1:ceil(length(data_in)/2)) = 1;
    bootstrap_ts = logical(bootstrap_ts(randperm(length(data_in))));
    bootstrap_ts(keep_idx == 0) = 0;
    
    
    % Compute the actual tuning curve using a bootstrapped sample
    if isempty(bins2) % 1d methods
        [actual_bs_MI(iShuff), ~, ~, actual_bs_prob_being_active(iShuff), actual_bs_TC(:,iShuff)] = MS_get_spatial_information(data_in(bootstrap_ts),this_shuff_samples(bootstrap_ts), pos_in(bootstrap_ts,1), bins1);
        
        % Compute the shuffled tuning curve using the same bootstrapped sample
        [shuffled_bs_MI(iShuff), ~, ~, shuffled_bs_prob_being_active(iShuff), shuffled_bs_TC(:,iShuff)]  = MS_get_spatial_information(data_in(bootstrap_ts),this_shuff_samples(bootstrap_ts), pos_in(bootstrap_ts,1), bins1);
        
    else % 2d methods.
        %         [actual_bootstrap_MI(iShuff), actual_bootstrap_PDF(:,iShuff), ~, actual_bootstrap_prob_being_active(iShuff), actual_bootstrap_tuning_curve(:,iShuff) ] = data_in_get_spatial_information(data_in.Binary(:,iC), behav_aligned.position(:,1), bins, split_ts);
        [actual_bs_MI(iShuff), ~, ~, actual_bs_prob_being_active(iShuff), actual_bs_TC(:,:,iShuff)] = MS_get_spatial_information_2D(data_in(bootstrap_ts), pos_in(bootstrap_ts,:), bins1,bins2);
        
        % Compute the shuffled tuning curve using the same bootstrapped sample
        [shuffled_bs_MI(iShuff),~, ~, shuffled_bs_prob_being_active(iShuff), shuffled_bs_TC(:,:,iShuff)] = MS_get_spatial_information_2D(this_shuff(bootstrap_ts),pos_in(bootstrap_ts,:), bins1,bins2);
    end
    
    
end % shuffle

%% get the stats for output
% Find the 95% confidence interval
sorted_BS_tuning_curves = sort(actual_bs_TC,2);
CI_idx_loc = 0.95*nShuff/2;
median_idx = round(nShuff/2);
upper_CI95_idx = median_idx+CI_idx_loc;
lower_CI95_idx = median_idx-CI_idx_loc;

% This will make sure that upper and lower bounds are withing the actual bootstrap data
upper_CI95_idx(upper_CI95_idx > nShuff) = nShuff;
upper_CI95_idx(upper_CI95_idx < 1) = 1;
lower_CI95_idx(lower_CI95_idx > nShuff) = nShuff;
lower_CI95_idx(lower_CI95_idx < 1) = 1;


%% gather for output.
stats_out.nShuff = nShuff;

if isempty(bins2)
    stats_out.mean = mean(actual_bs_TC,2);
    stats_out.median = median(actual_bs_TC,2);
else
    stats_out.mean = mean(actual_bs_TC, 3);
    stats_out.median = median(actual_bs_TC, 3);
end

stats_out.upper_CI95 = sorted_BS_tuning_curves(:,upper_CI95_idx);
stats_out.lower_CI95 = sorted_BS_tuning_curves(:,lower_CI95_idx);

stats_out.actual_TC = actual_bs_TC;
stats_out.actual_MI = actual_bs_MI;
stats_out.actual_p_active = actual_bs_prob_being_active;

stats_out.shuff_TC = shuffled_bs_TC;
stats_out.shuff_MI = shuffled_bs_MI;
stats_out.shuff_p_active = shuffled_bs_prob_being_active;




end % function