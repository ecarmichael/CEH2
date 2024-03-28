function [behav_enc, behav_rec] = CE_DLC_Rad(data_dir, inter_dir)
%% CE_DLC_Rad: collects the DLC outputs from the encoding and recall phases of the Radial task







%% check for behav recordings in the main dir


recs = dir([data_dir filesep '*_*']);

% load all the behav files
behav = [];
for iR = length(recs):-1:1
    
    cd([recs(iR).folder filesep recs(iR).name])
    
    if exist([recs(iR).folder filesep recs(iR).name filesep 'minicam1'], 'dir')
        
        [~, behav{iR}] = MS_DLC2TSD([recs(iR).folder filesep recs(iR).name filesep 'minicam1'], [], [9.33 9.33]);
        keep_idx(iR) = 1;
        rec_len(iR) = behav{iR}.time(end)  -  behav{iR}.time(1);
        rec_time(iR) = behav{iR}.json.recordingStartTime.msecSinceEpoch;
        
    else
        keep_idx(iR) = 0;
        rec_len(iR) = 0;
        rec_time(iR) = 0;
    end
    
end

behav(~keep_idx) = [];
rec_len(~keep_idx) = [];
rec_time(~keep_idx) = [];

if length(rec_len) > 2
    min_len = 600;
    
    rm_idx = rec_len < min_len;
    
    behav(rm_idx) = [];
    rec_len(rm_idx) = [];
    rec_time(rm_idx) = [];
    
end

    [~, enc_idx] = min(rec_time);
    [~, rec_idx] = max(rec_time);

behav_enc = behav{enc_idx}; 
behav_rec = behav{rec_idx}; 

fprintf('Encoding phase started at %2.0f:%2.0f:%2.0f (dur = %0.2fmins).\n', behav_enc.json.recordingStartTime.hour,behav_enc.json.recordingStartTime.minute, behav_enc.json.recordingStartTime.second, (behav_enc.time(end) - behav_enc.time(1))/60)
fprintf('Recall phase started at %2.0f:%2.0f:%2.0f (dur = %0.2fmins).\n', behav_rec.json.recordingStartTime.hour,behav_rec.json.recordingStartTime.minute, behav_rec.json.recordingStartTime.second, (behav_rec.time(end) - behav_rec.time(1))/60)



save([inter_dir filesep 'behav_enc.mat'], 'behav_enc', '-v7.3')
save([inter_dir filesep 'behav_rec.mat'], 'behav_rec', '-v7.3')

end

