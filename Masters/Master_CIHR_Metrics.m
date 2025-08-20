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


%% genotype and virus arrays



%% NOL %%%%%%%%

data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\CIHR_NOL_2025'; 
%%
cd(data_dir)

% load the opto_grant.xlsx file containing the IDs, genotypes, and viruses. 

opto_tbl = readtable("Opto_grant.xlsx", 'VariableNamingRule','preserve'); 


% get all of the "*markers.csv' files for the scored interactions

mkr_list = dir('*markers.csv'); 

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


    tbl = readtable(mkr_list(ii).name); 

    % get all the object epochs
    obj_1_idx = find((tbl.object_id == 1 ) & contains(tbl.marker_type, 'start')); 
    obj_2_idx = find((tbl.object_id == 2 ) & contains(tbl.marker_type, 'start')); 

    % remove any blocks starting after 600seconds
    rm_idx = tbl.marker_time(obj_1_idx) > 600; 
    obj_1_idx(rm_idx) = []; 
    
    rm_idx = tbl.marker_time(obj_2_idx) > 600; 
    obj_2_idx(rm_idx) = []; 

    % loop over events and get the duration
    obj_1 = []; 
    for jj = length(obj_1_idx):-1:1

        if  tbl.marker_time(obj_1_idx(jj)+1) > 600
            obj_1(jj) = 600 - tbl.marker_time(obj_1_idx(jj));
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
    data_out{ii, 5} = length(obj_1); 
    data_out{ii, 6} = sum(obj_1); 
    data_out{ii, 7} = length(obj_2); 
    data_out{ii, 8} = sum(obj_2); 
    data_out{ii, 9} = (sum(obj_2) - sum(obj_1))/(sum(obj_2)+sum(obj_1)); 
    data_out{ii, 10} = (length(obj_2) - length(obj_1))/(length(obj_2)+length(obj_1)); 

    fprintf("%s - %s Obj 1 n: %.0f  t:%.2f    |  Obj 2 n: %.0f  t:%.2f   | <strong>%.1f</strong> <strong>(%.1f)</strong>\n", Sub{ii}, Sess{ii}, length(obj_1), sum(obj_1), length(obj_2), sum(obj_2), data_out{ii,7}, data_out{ii,8})


end

% convert to table

tbl_out = cell2table(data_out, "VariableNames",{'Subject', 'Session','Geno', 'Virus', 'Obj1_n', 'Obj1_t', 'Obj2_n', 'Obj2_t', 'DI_n', 'DI_t'});


%% plot
c_ord = MS_linspecer(8); 

e_idx = contains(tbl_out.Session, 'ENC');

r_idx = contains(tbl_out.Session, 'REC'); 

ctrl_idx =  contains(tbl_out.Virus, 'Sham'); 

cheta_idx =  contains(tbl_out.Virus, 'Cheta'); 

bi_idx =  contains(tbl_out.Virus, 'Bipole');

art_idx =  contains(tbl_out.Virus, 'ArchT');

sst_idx =  contains(tbl_out.Geno, 'Sst');

figure(191)
clf
% overall
subplot(2,2,1:2)
[hb, eb, sc, p, stats] = MS_bar_w_err(tbl_out.DI_t(e_idx)', tbl_out.DI_t(r_idx)', [c_ord(1,:);c_ord(end,:)],1,  'ttest2', [1 2]);

% 
[hb, eb, sc, p, stats] = MS_bar_w_err(tbl_out.DI_t(e_idx & ctrl_idx)', tbl_out.DI_t(r_idx & ctrl_idx)', [c_ord(3,:);c_ord(3,:)],1,  'ttest2', [4 5]);

[hb, eb, sc, p, stats] = MS_bar_w_err(tbl_out.DI_t(e_idx &  bi_idx)', tbl_out.DI_t(r_idx &  bi_idx)', [c_ord(4,:);c_ord(4,:)],1,  'ttest2', [7 8]);

[hb, eb, sc, p, stats] = MS_bar_w_err(tbl_out.DI_t(e_idx & (art_idx))', tbl_out.DI_t(r_idx & (art_idx))', [c_ord(7,:);c_ord(8,:)],1,  'ttest2', [10 11]);

[hb, eb, sc, p, stats] = MS_bar_w_err(tbl_out.DI_t(e_idx & cheta_idx)', tbl_out.DI_t(r_idx & cheta_idx)', [c_ord(5,:);c_ord(5,:)],1,  'ttest2', [13 14]);

set(gca, 'XTick', [1 2 4 5 7 8 10 11 13 14], 'XTickLabel', {'All - Enc', 'All - Rec', 'Ctrl - Enc', 'Ctrl - Rec', 'Bipole - Enc', 'Bipole - Rec'...
    'ArchT - Enc', 'ArchT - Rec', 'Cheta - Enc', 'Cheta - Rec'}, 'XTickLabelRotation', 45)
subplot(2,2,3)
[hb, eb, sc, p, stats] = MS_bar_w_err(tbl_out.DI_t(r_idx & ctrl_idx)', tbl_out.DI_t(r_idx & (art_idx | bi_idx))', [[ .6 .6 .6];c_ord(1,:)],1,  'ttest2', [1 2]);

hb.FaceColor = "none"; 
sc{1}.MarkerFaceColor = [ .6 .6 .6]; 
sc{2}.MarkerFaceColor = c_ord{1}; 


%%%%%% to do %%%
%convert into vectors for paired tests. 

% plot as all four conditions for encoding  and then for recall

% plot as all four conditions and comparisons

% make a gittered scatter with the mouse labels to see if any odd sessions
% are outliers. 

%% get the distance traveled and occupancy heat maps

data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\CIHR_NOL_2025'; 
cd(data_dir)
% get all of the "*markers.csv' files for the scored interactions

dlc_list = dir('*shuffle*.csv'); 

% loop over and collect the interactions 
Sub = cell(length(dlc_list),1); Sess = Sub; 
data_out = cell(length(dlc_list), 6); 

dist = []; c_dist = {}; occ_mat = []; 

for ii  = length(dlc_list):-1:1

    vname = [dlc_list(ii).name(1:strfind(dlc_list(ii).name, 'DLC')-1) '.mp4'];

    pos = MS_DLC2TSD_single(dlc_list(ii).name, vname, [805/38 805/38]); 

    % get an occupancy map
    [~, occ_mat(:,:,ii)] = MS_event_rate_map(zeros(1,size(pos.tvec,2)), pos.tvec, pos.data(1:2,:)', 2.5,0, 2, 24:1:68,4:1:46 ); 

    % get the total distance traveled. 
    dist(ii) = trapz(pos.tvec, pos.data(end-1,:)); 
    cdist{ii} = cumtrapz(pos.tvec, pos.data(end-1,:)); 
end