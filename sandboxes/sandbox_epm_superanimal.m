%% sandbox_EPM_superanimal

fnames = dir('*.csv');
vids = dir('*.mp4');


f_list = []; v_list= []; 
for ii = 1:length(fnames)

    f_list{ii} = fnames(ii).name; 
    v_list{ii} = vids(ii).name; 

    % vid_name = contains(fname)
    % [pos] = MS_DLC2TSD_single(fnames(ii).name)
    
end


%% run the DLC

for ii = 1:length(f_list)

    splt_idx = strfind(f_list{ii}, '_super'); 

    vid_name = [f_list{ii}(1:splt_idx-1) '.mp4'];

    [pos] = MS_DLC2TSD_single(f_list{ii}, vid_name)

end