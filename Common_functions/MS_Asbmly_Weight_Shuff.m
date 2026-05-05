function [shuff_stats, w_thresh] = MS_Asbmly_Weight_Shuff(P_temp, data_h, tvec, nShuff)
%% MS_Asbmbly_Weight_Shuff: apply a weight permutation against the real data. 

if nargin < 4
    nShuff = 100; % Default value for nShuff if not provided
end

%% run the shuffle

tic % time counter
% loop over assembly templates
for iA = size(P_temp,2):-1:1
    for iS = nShuff:-1:1
        % Shuffle the data
        S_temp(:,iA,iS) = P_temp(randperm(size(P_temp, 1)), iA);

        S_prog = assembly_activity(S_temp(:,iA,iS),data_h');

        wake_shuff_mat(:,iA,iS) =  S_prog(1,:);

    end
w_thresh(iA) = prctile(wake_shuff_mat(:,iA,:), 99, 'all');

n_sig_R(iA,iS) = sum(wake_shuff_mat(:,iA,:) > w_thresh(iA)); 
end

% collect 
t_dur = (tvec(end) - tvec(1)) / 60; 

shuff_stats.shuff_proj = wake_shuff_mat; 
shuff_stats.w_thresh  = w_thresh; 
shuff_stats.w_thresh_all  = prctile(wake_shuff_mat, 99, 'all'); 


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


% print it
fprintf('Weight Shuff mean R thresh: %.2f [elapsed time: %.2fs]\n', shuff_stats.w_thresh_all, toc)