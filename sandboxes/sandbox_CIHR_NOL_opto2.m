%% sandbox_CIHR_NOL_2

% gathers the NOL and TFC metrics and makes some nice plots. 

% NOL up top | TFC below 
%[Note Obj 1 is constant and Obj 2 moves 
%       _____________
%       |           |
%       |        2* |
%       X           |
%       |  1     2  |
%       |___________|

%%%%%%%%%%% chronate M7-m9 recall

%% NOL %%%%%%%%
if ispc
data_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\CIHR_NOL_2025\Opto_cohort2'; 
else
data_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/CIHR_2025/CIHR_NOL_2025'; 
end
%
cd(data_dir)

% load the opto_grant.xlsx file containing the IDs, genotypes, and viruses. 

opto_tbl = readtable("Opto_grant.xlsx", 'VariableNamingRule','preserve'); 


% get all of the "*markers.csv' files for the scored interactions

mkr_list = dir('*markers.csv'); 

for ii = 1:length(mkr_list)
    if (sum(contains(mkr_list(ii).name, 'p1')) > 0) || ( sum(contains(mkr_list(ii).name, 'p2')) > 0)
        rm_idx(ii) = true;
    else
        rm_idx(ii) = false;
    end
end

mkr_list(rm_idx) = []; 

% loop over and collect the interactions 
Sub = cell(length(mkr_list),1); Sess = Sub; Geno = Sub; Virus = Sub; 

data_out = cell(length(mkr_list), 6); 
for ii  = 1:length(mkr_list)
    
    s_idx = strfind(mkr_list(ii).name, '_'); 
    Sub{ii} = mkr_list(ii).name(1:s_idx(1)-1); 
    Sess{ii} = upper(mkr_list(ii).name(s_idx(end-1)+1:s_idx(end)-1));

    % pull the genotype and virus
    this_idx = find(contains(opto_tbl.Subject, Sub{ii})); 
    Geno{ii} = opto_tbl.Strain{this_idx}; 
    Virus{ii} = opto_tbl.Injection{this_idx}; 
    Sex{ii} = opto_tbl.Sex{this_idx}; 

    tbl = readtable(mkr_list(ii).name); 

    % get all the object epochs
    obj_1_idx = find((tbl.object_id == 1 ) & contains(tbl.marker_type, 'start')); 
    obj_2_idx = find((tbl.object_id == 2 ) & contains(tbl.marker_type, 'start')); 

    % remove any blocks starting after 600seconds
    rm_idx = tbl.marker_time(obj_1_idx) > 300; 
    obj_1_idx(rm_idx) = []; 
    
    rm_idx = tbl.marker_time(obj_2_idx) > 300; 
    obj_2_idx(rm_idx) = []; 

    % loop over events and get the duration
    obj_1 = []; 
    for jj = length(obj_1_idx):-1:1

        if  tbl.marker_time(obj_1_idx(jj)+1) > 300
            obj_1(jj) = 300 - tbl.marker_time(obj_1_idx(jj));
        else
            obj_1(jj) = tbl.marker_time(obj_1_idx(jj)+1) - tbl.marker_time(obj_1_idx(jj)); 
        end
    end

    obj_2 = []; 
    for jj = length(obj_2_idx):-1:1
        obj_2(jj) = tbl.marker_time(obj_2_idx(jj)+1) - tbl.marker_time(obj_2_idx(jj)); 
    end

    % collect in a sheet
    data_out{ii, 1} = Sub{ii}; 
    data_out{ii, 2} = Sess{ii}; 
    data_out{ii, 3} = Geno{ii}; 
    data_out{ii, 4} = Virus{ii}; 
    data_out{ii, 5} = Sex{ii};
    data_out{ii, 6} = length(obj_1); 
    data_out{ii, 7} = sum(obj_1); 
    data_out{ii, 8} = length(obj_2); 
    data_out{ii, 9} = sum(obj_2); 
    data_out{ii, 10} = (sum(obj_2) - sum(obj_1))/(sum(obj_2)+sum(obj_1)); 
    data_out{ii, 11} = (length(obj_2) - length(obj_1))/(length(obj_2)+length(obj_1)); 

    fprintf("%s - %s Obj 1 n: %.0f  t:%.2f    |  Obj 2 n: %.0f  t:%.2f   | <strong>%.1f</strong> <strong>(%.1f)</strong>\n", Sub{ii}, Sess{ii}, length(obj_1), sum(obj_1), length(obj_2), sum(obj_2), data_out{ii,10}, data_out{ii,11})


end



% convert to table

tbl_out = cell2table(data_out, "VariableNames",{'Subject', 'Session','Geno', 'Virus','Sex', 'Obj1_n', 'Obj1_t', 'Obj2_n', 'Obj2_t', 'DI_n', 'DI_t'});

% make a paired table as well. 

enc_idx = (contains(data_out(:,2), 'ENC')); 

data_pair_e = data_out(:,[1 2 3 4 5 6 7 8 9 10 11]); 
data_pair_r = data_out(:,[1 2 3 4 5 6 7 8 9 10 11]); 
data_pair_e(~enc_idx,:) = []; 
data_pair_r(enc_idx,:) = []; 

data_pair = [data_pair_e, data_pair_r(:, [1 2 6 7 8 9 10 11])];
%%%%%% double check this is getting the right variables %%%%%%%%%%%%
tbl_pairs = cell2table(data_pair, "VariableNames",{'Subject', 'Session','Geno', 'Virus','Sex','E_Obj1_t', 'E_Obj1_n', 'E_Obj2_t', 'E_Obj2_n' 'E_DI_n', 'E_DI_t','SubRr', 'SessR','R_Obj1_t', 'R_Obj1_n', 'R_Obj2_t', 'R_Obj2_n', 'R_DI_n', 'R_DI_t'});

% check the subjects line up. 
for ii = 1:size(tbl_pairs.Subject,1)
    if ~strcmp(tbl_pairs.Subject, tbl_pairs.SubRr)
        disp([tbl_pairs.Subject(ii) ' - ' tbl_pairs.SubRr(ii)]);
    end

    if strcmp(tbl_pairs.Session, tbl_pairs.SessR)
        disp([tbl_pairs.Session(ii) ' - ' tbl_pairs.SessR(ii)]);
    end

    if (sum([tbl_pairs.E_Obj1_t(ii),tbl_pairs.E_Obj2_t(ii)]) < 5) || (sum([tbl_pairs.R_Obj1_t(ii),tbl_pairs.R_Obj2_t(ii)]) < 5)
        fprintf('%s  E = %.2fs  | R = %.2fs \n', tbl_pairs.Subject{ii}, sum([tbl_pairs.E_Obj1_t(ii),tbl_pairs.E_Obj2_t(ii)]),sum([tbl_pairs.R_Obj1_t(ii),tbl_pairs.R_Obj2_t(ii)])); 
    end
end
