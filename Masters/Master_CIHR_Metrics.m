%% Master_CIHR_Metrics

% gathers the NOL and TFC metrics and makes some nice plots. 

% NOL up top | TFC below 
%[Note Obj 1 is constant and Obj 2 moves 
%       _____________
%       |           |
%       |        2* |
%       X           |
%       |  1     2  |
%       |___________|




%% NOL %%%%%%%%

data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\CIHR_NOL_2025'; 
cd(data_dir)
% get all of the "*markers.csv' files for the scored interactions

mkr_list = dir('*markers.csv'); 

% loop over and collect the interactions 
Sub = cell(length(mkr_list),1); Sess = Sub; 
data_out = cell(length(mkr_list), 6); 
for ii  = 1:length(mkr_list)
    
    s_idx = strfind(mkr_list(ii).name, '_'); 
    Sub{ii} = mkr_list(ii).name(1:s_idx(1)-1); 
    Sess{ii} = mkr_list(ii).name(s_idx(end-1)+1:s_idx(end)-1);

    tbl = readtable(mkr_list(ii).name); 

    % get all the object epochs
    obj_1_idx = find((tbl.object_id == 1 ) & contains(tbl.marker_type, 'start')); 
    obj_2_idx = find((tbl.object_id == 2 ) & contains(tbl.marker_type, 'start')); 

    % loop over events and get the duration
    obj_1 = []; 
    for jj = length(obj_1_idx):-1:1
        obj_1(jj) = tbl.marker_time(obj_1_idx(jj)+1) - tbl.marker_time(obj_1_idx(jj)); 
    end

    obj_2 = []; 
    for jj = length(obj_2_idx):-1:1
        obj_2(jj) = tbl.marker_time(obj_2_idx(jj)+1) - tbl.marker_time(obj_2_idx(jj)); 
    end

    fprintf("%s - %s Obj 1 n: %.0f  t:%.2f    |  Obj 2 n: %.0f  t:%.2f\n", Sub{ii}, Sess{ii}, length(obj_1), sum(obj_1), length(obj_2), sum(obj_2))
    % collect in a sheet
    data_out{ii, 1} = Sub{ii}; 
    data_out{ii, 2} = Sess{ii}; 
    data_out{ii, 3} = length(obj_1); 
    data_out{ii, 4} = sum(obj_1); 
    data_out{ii, 5} = length(obj_2); 
    data_out{ii, 6} = sum(obj_2); 

end

%% get the distance traveled and occupancy heat maps

data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\CIHR_NOL_2025'; 
cd(data_dir)
% get all of the "*markers.csv' files for the scored interactions

dlc_list = dir('*shuffle*.csv'); 

% loop over and collect the interactions 
Sub = cell(length(dlc_list),1); Sess = Sub; 
data_out = cell(length(dlc_list), 6); 

for ii  = 1:length(dlc_list)

    pos = MS_DLC2TSD_single(dlc_list(ii).name, )

end