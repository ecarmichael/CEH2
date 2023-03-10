%% sandbox_replay co-activity


% get data
if isunix
    load('/home/williamslab/Desktop/pv1060/HATD5/all_binary_post_REM.mat');
    load('/home/williamslab/Dropbox (Williams Lab)/Decoding/pv1060/HATD5/1000shuffling.mat'); 
end

clearvars -except all_binary_post_REM Final_start_frame Final_Replay_score

%% remove overlapping replays

[~, keep_idx] = findpeaks(diff(Final_start_frame), 'MinPeakHeight', 5);
keep_idx = [1, keep_idx(1:end)-1]; 
%% get the replay indicies and make some PETHs

win_s = 33*2; 

for ii = 1:length(Final_start_frame)

%     figure(ii)
    this_data = all_binary_post_REM(Final_start_frame(ii) - win_s:Final_start_frame(ii) + win_s,:);
%     MS_Ca_Raster(this_data', Final_start_frame(ii) - win_s:Final_start_frame(ii) + win_s)
    
    data_mat(:,:,ii) = this_data; 
    
    
end


%% plot the Peri-replay activity matrix

% imagesc(