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
w_thresh(iA) = prctile(squeeze(wake_shuff_mat(:,iA,:)), 99, 'all');

n_sig_R(iA,:) = sum(squeeze(wake_shuff_mat(:,iA,:)) > w_thresh(iA),1); 
end

% collect 
t_dur = (tvec(end) - tvec(1)) / 60; 

shuff_stats.shuff_proj = wake_shuff_mat; 
shuff_stats.w_thresh  = w_thresh; 
shuff_stats.w_thresh_all  = prctile(wake_shuff_mat, 99, 'all'); 

% print it
fprintf('Weight Shuff mean R thresh: %.2f [elapsed time: %.2fs]\n', shuff_stats.w_thresh_all, toc)