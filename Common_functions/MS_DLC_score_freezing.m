function  out = MS_DLC_score_freezing(fname)
%% MS_score_freezing: score freezing based on movement in DLC tracking data. if multiple body parts, average across to get a better measure. 









%% initialize 

% find the label file
f_list = dir('*filtered.csv');

keep_idx = zeros(1, length(f_list)); 

for f = 1:length(f_list)
    if contains(f_list(f).name, fname(1:end-4))
        keep_idx(f) = true; 
    else
        keep_idx(f) = false;
    end
end

keep_idx = logical(keep_idx); 

dlc_name = f_list(keep_idx).name; 
dlc_dir = f_list(keep_idx)

pos = MS_DLC2TSD_single(dlc_name, fname); 