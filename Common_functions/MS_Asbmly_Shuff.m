function [shuff_stats, w_thresh] = MS_Asbmly_Shuff(data_h, tvec, nShuff, opts)
%% MS_Asbmbly_Shuff: apply the shuffling checks only.  Used to recompute the shuffle distributions post-hoc (eg: if a large shuffle number is needed). 

if nargin < 3
    nShuff = 100; % Default value for nShuff if not provided
        opts.threshold.method = 'MarcenkoPastur';
    opts.Patterns.method = 'ICA';
    opts.Patterns.number_of_iterations = 500;
    opts.threshold.permutations_percentile= 95;
    opts.threshold.number_of_permutations= 500;
elseif nargin < 4
    opts.threshold.method = 'MarcenkoPastur';
    opts.Patterns.method = 'ICA';
    opts.Patterns.number_of_iterations = 500;
    opts.threshold.permutations_percentile= 95;
    opts.threshold.number_of_permutations= 500;
end

%% run the shuffle
tic
% places for data to go
wake_shuff_mat = [];
all_shuff = []; 

Ass_shuff = NaN(1,nShuff);
for iS = nShuff:-1:1
    tic
    shuff_data = NaN(size(data_h));
    for ic = 1:size(data_h,2)
        shuff_data(:,ic) = circshift(data_h(:,ic), floor(MS_randn_range(1,1,1,size(data_h,1))));
    end
    
    this_ass = assembly_patterns(shuff_data', opts);
    if ~isempty(this_ass)
        S_prog = assembly_activity(this_ass,shuff_data');
        
        wake_shuff_mat(iS,:) =  S_prog(1,:);
        keep_idx(iS) = 1;
    else
        wake_shuff_mat(iS,:) = NaN(1,length(shuff_data));
        keep_idx(iS) = 0;
    end
    
    if sum(max(this_ass) > 0.2) >0
        Ass_shuff(iS) = sum(max(this_ass) > 0.2);
    else
        Ass_shuff(iS) = 0;
    end
    all_shuff(:,:,iS) = shuff_data; 
end


w_thresh = prctile(wake_shuff_mat(wake_shuff_mat >0), 99, 'all');
t_dur = (tvec(end) - tvec(1)) / 60; 

shuff_stats.shuff_data = all_shuff; 
shuff_stats.shuff_proj = wake_shuff_mat; 
shuff_stats.shuff_n = Ass_shuff; 
shuff_stats.mean = mean(Ass_shuff); 
shuff_stats.sd = std(Ass_shuff); 
shuff_stats.shuff_r = Ass_shuff./t_dur; 
shuff_stats.mean_r = mean(Ass_shuff./t_dur); 
shuff_stats.sd_r = std(Ass_shuff./t_dur); 
shuff_stats.p95 = prctile(Ass_shuff, 95, 'all');
shuff_stats.p99 = prctile(Ass_shuff, 99, 'all');
shuff_stats.p95_r = prctile(Ass_shuff./t_dur, 95, 'all');
shuff_stats.p99_r = prctile(Ass_shuff./t_dur, 99, 'all');
shuff_stats.w_thresh  = w_thresh; 

fprintf('Shuff mean assemblies %.2f +/-%.2f | Rate:  %.2f +/-%.2f/min  | R thresh: %.2f [elapsed time: %.2fs]\n', shuff_stats.mean, shuff_stats.sd, shuff_stats.mean_r, shuff_stats.sd_r, w_thresh, toc)


