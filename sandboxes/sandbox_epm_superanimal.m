%% sandbox_EPM_superanimal

fnames = dir('*.csv');
vids = dir('*.mp4');


dir_list = []; f_list = []; v_list= []; 
for ii = 1:length(fnames)
    
    dir_list{ii} = fnames(ii).folder; 
    f_list{ii} = fnames(ii).name; 
    v_list{ii} = vids(ii).name; 

    % vid_name = contains(fname)
    % [pos] = MS_DLC2TSD_single(fnames(ii).name)
    
end


%% run the DLC
% boxes = []; 

for ii = 7:length(f_list)

    splt_idx = strfind(f_list{ii}, '_super');

    vid_name = [f_list{ii}(1:splt_idx-1) '.mp4'];

    if isempty(boxes)
        [emp_idx, labels, boxes] = MS_DLC_EPM_super(f_list{ii},'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\PoxR1\EPM\inter_temp' , [], 0);
    else
        [emp_idx, labels] = MS_DLC_EPM_super(f_list{ii},'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\PoxR1\EPM\inter_temp' , boxes, 0);
    end
    % [pos] = MS_DLC2TSD_single(f_list{ii}, vid_name);

end