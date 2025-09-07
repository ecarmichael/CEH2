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

%%%%%%%%%%% chronate M7-m9 recall

%% NOL %%%%%%%%

data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\CIHR_NOL_2025'; 
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

data_pair = [data_pair_e, data_pair_r(:, [1 2 6 7 8 9])];
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




%% reshash mice with lost optics
f_idx  = find(contains(tbl_out.Subject, 'M9'));

for ii = 1:length(f_idx)
    tbl_out.Virus{f_idx(ii)} = 'Sham';
end

f_idx  = contains(tbl_out.Subject, {'M23' 'M4'});
tbl_out(f_idx,:) = []; 


f_idx  = find(contains(tbl_pairs.Subject, 'M9'));

for ii = 1:length(f_idx)
    tbl_pairs.Virus{f_idx(ii)} = 'Sham';
end

f_idx  = contains(tbl_pairs.Subject, {'M23' 'M4'});
tbl_pairs(f_idx,:) = []; 

%% plot
c_ord = MS_linspecer(8); 

ctrl_idx =  contains(tbl_pairs.Virus, 'Sham'); 

cheta_idx =  contains(tbl_pairs.Virus, 'Cheta'); 

bi_idx =  contains(tbl_pairs.Virus, 'Bipole');

art_idx =  contains(tbl_pairs.Virus, 'ArchT');

sst_idx =  contains(tbl_pairs.Geno, 'Sst');

figure(192)
this_data_e = tbl_pairs.E_DI_t; 
this_data_r = tbl_pairs.R_DI_t; 

clf
% overall
subplot(2,2,1)
[hb, eb, sc, p, stats] = MS_bar_w_err(this_data_e', this_data_r', [c_ord(1,:);c_ord(end,:)],1,  'ttest', [1 2]);

% 
[hb, eb, sc, p, stats] = MS_bar_w_err(this_data_e(ctrl_idx)', this_data_r(ctrl_idx)', [.6 .6 .6 ;.6 .6 .6],1,  'ttest', [4 5]);


[hb, eb, sc, p, stats] = MS_bar_w_err(this_data_e(bi_idx)', this_data_r(bi_idx)', [c_ord(2,:);c_ord(2,:)],1,  'ttest', [7 8]);


[hb, eb, sc, p, stats] = MS_bar_w_err(this_data_e(art_idx)', this_data_r(art_idx)', [c_ord(3,:);c_ord(3,:)],1,  'ttest', [10 11]);

[hb, eb, sc, p, stats] = MS_bar_w_err(this_data_e(cheta_idx)', this_data_r(cheta_idx)', [c_ord(5,:);c_ord(5,:)],1,  'ttest', [13 14]);

set(gca, 'XTick', [1 2 4 5 7 8 10 11 13 14], 'XTickLabel', {'All - Enc', 'All - Rec', 'Ctrl - Enc', 'Ctrl - Rec', 'Bipole - Enc', 'Bipole - Rec'...
    'ArchT - Enc', 'ArchT - Rec', 'Cheta - Enc', 'Cheta - Rec'}, 'XTickLabelRotation', 45)

ylabel({'Discrimination index'; '(Obj2 - Obj1)/(Obj1+obj2)'})

% subplot(2,2,3)
% [hb, eb, sc, p, stats] = MS_bar_w_err(this_data(r_idx & ctrl_idx)', this_data(r_idx & (art_idx | bi_idx))', [[ .6 .6 .6];c_ord(1,:)],1,  'ttest2', [1 2]);


% plot the interaction number and events per type as a scatter
subplot(2,2,2)
hold on
scatter(tbl_pairs.R_DI_t(ctrl_idx), tbl_pairs.R_DI_n(ctrl_idx), 100, 'MarkerFaceColor', [.6 .6 .6],"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
scatter(tbl_pairs.R_DI_t(bi_idx), tbl_pairs.R_DI_n(bi_idx), 100, 'MarkerFaceColor', c_ord(2,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
scatter(tbl_pairs.R_DI_t(art_idx), tbl_pairs.R_DI_n(art_idx), 100, 'MarkerFaceColor', c_ord(3,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
scatter(tbl_pairs.R_DI_t(cheta_idx), tbl_pairs.R_DI_n(cheta_idx), 100, 'MarkerFaceColor', c_ord(5,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
xline(0); yline(0);
ylabel('DI: n interactions')
xlabel('DI: time')
legend({'Control', 'Bipole', 'ArchT', 'Cheta'}, 'Box','off', 'Location','northwest')
xlim([-1 1]); ylim([-1 1])

subplot(2,2,4)
hold on
scatter(tbl_pairs.R_Obj1_n(ctrl_idx), tbl_pairs.R_Obj2_n(ctrl_idx), 100, 'MarkerFaceColor', [.6 .6 .6],"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
scatter(tbl_pairs.R_Obj1_n(bi_idx), tbl_pairs.R_Obj2_n(bi_idx), 100, 'MarkerFaceColor', c_ord(2,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
scatter(tbl_pairs.R_Obj1_n(art_idx), tbl_pairs.R_Obj2_n(art_idx), 100, 'MarkerFaceColor', c_ord(3,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
scatter(tbl_pairs.R_Obj1_n(cheta_idx), tbl_pairs.R_Obj2_n(cheta_idx), 100, 'MarkerFaceColor', c_ord(5,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
xline(0); yline(0);
ylabel('DI: n interactions')
xlabel('DI: time')
legend({'Control', 'Bipole', 'ArchT', 'Cheta'}, 'Box','off', 'Location','northwest')
xlim([-1 1]); ylim([-1 1])

subplot(2,2,3)
cla
hold on
yline(0)
this_data = tbl_pairs.R_DI_t;

% ctrl
idx = find(ctrl_idx);
x_vec = sort(MS_randn_range(1,sum(ctrl_idx), .8, 1.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', [.6 .6 .6],"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
for ii = 1:length(idx) 
    if contains(tbl_pairs.Sex(idx(ii)), 'M')
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)))
    end
    disp(tbl_pairs.Subject(idx(ii)))
end


% bipole
idx = find(bi_idx & ~sst_idx);
x_vec = sort(MS_randn_range(1,sum( bi_idx & ~sst_idx), 1.8, 2.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', c_ord(2,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
for ii = 1:length(idx) 
    if contains(tbl_pairs.Sex(idx(ii)), 'M')
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)))
    end
end

idx = find(bi_idx & sst_idx);
x_vec = sort(MS_randn_range(1,sum( bi_idx & sst_idx), 1.8, 2.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', c_ord(2,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"diamond")
for ii = 1:length(idx) 
    if contains(tbl_pairs.Sex(idx(ii)), 'M')
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)))
    end
end

% archt
idx = find(art_idx & ~sst_idx);
x_vec = sort(MS_randn_range(1,sum(art_idx& ~sst_idx), 2.8, 3.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', c_ord(3,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
for ii = 1:length(idx) 
    if contains(tbl_pairs.Sex(idx(ii)), 'M')
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)))
    end
end

idx = find(art_idx & sst_idx);
x_vec = sort(MS_randn_range(1,sum( art_idx& sst_idx), 2.8, 3.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', c_ord(3,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"diamond")
for ii = 1:length(idx) 
    if contains(tbl_pairs.Sex(idx(ii)), 'M')
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)))
    end
end

idx = find(cheta_idx & ~sst_idx);
x_vec = sort(MS_randn_range(1,sum( cheta_idx & ~sst_idx),3.8, 4.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', c_ord(5,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
for ii = 1:length(idx)
    if contains(tbl_pairs.Sex(idx(ii)), 'M')
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)))
    end
end

idx = find( cheta_idx & sst_idx);
x_vec = sort(MS_randn_range(1,sum(cheta_idx & sst_idx),3.8, 4.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', c_ord(5,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"diamond")
for ii = 1:length(idx)
    if contains(tbl_pairs.Sex(idx(ii)), 'Male')
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_pairs.Subject(idx(ii)))
    end
end

set(gca, 'XTick', 1:4, 'XTickLabel', { 'Ctrl - Rec' 'Bipole - Rec'...
    'ArchT - Rec' 'Cheta - Rec'}, 'XTickLabelRotation', 45)

ylim([-1 1]); 
xlim ([-.5 4.5])
ylabel({'Discrimination index'; '(Obj2 - Obj1)/(Obj1+obj2)'})
%%%%%% to do %%%
%convert into vectors for paired tests. 

% plot as all four conditions for encoding  and then for recall

% plot as all four conditions and comparisons

% make a gittered scatter with the mouse labels to see if any odd sessions
% are outliers. 


%% plot draft without paired values.
c_ord = MS_linspecer(8); 

e_idx = contains(tbl_out.Session, 'ENC');

r_idx = contains(tbl_out.Session, 'REC'); 

ctrl_idx =  contains(tbl_out.Virus, 'Sham'); 

cheta_idx =  contains(tbl_out.Virus, 'Cheta'); 

bi_idx =  contains(tbl_out.Virus, 'Bipole');

art_idx =  contains(tbl_out.Virus, 'ArchT');

sst_idx =  contains(tbl_out.Geno, 'Sst');

figure(192)
this_data = tbl_out.DI_t; 

clf
% overall
subplot(2,2,1)
[hb, eb, sc, p, stats] = MS_bar_w_err(this_data(e_idx)', this_data(r_idx)', [c_ord(1,:);c_ord(end,:)],1,  'ttest2', [1 2]);
fprintf("")
% 
[hb, eb, sc, p, stats] = MS_bar_w_err(this_data(e_idx & ctrl_idx)', this_data(r_idx & ctrl_idx)', [.6 .6 .6 ;.6 .6 .6],1,  'ttest2', [4 5]);


[hb, eb, sc, p, stats] = MS_bar_w_err(this_data(e_idx &  bi_idx)', this_data(r_idx &  bi_idx)', [c_ord(2,:);c_ord(2,:)],1,  'ttest2', [7 8]);


[hb, eb, sc, p, stats] = MS_bar_w_err(this_data(e_idx & (art_idx))', this_data(r_idx & (art_idx))', [c_ord(3,:);c_ord(3,:)],1,  'ttest2', [10 11]);

[hb, eb, sc, p, stats] = MS_bar_w_err(this_data(e_idx & cheta_idx)', this_data(r_idx & cheta_idx)', [c_ord(5,:);c_ord(5,:)],1,  'ttest2', [13 14]);

set(gca, 'XTick', [1 2 4 5 7 8 10 11 13 14], 'XTickLabel', {'All - Enc', 'All - Rec', 'Ctrl - Enc', 'Ctrl - Rec', 'Bipole - Enc', 'Bipole - Rec'...
    'ArchT - Enc', 'ArchT - Rec', 'Cheta - Enc', 'Cheta - Rec'}, 'XTickLabelRotation', 45)

ylabel({'Discrimination index'; '(Obj2 - Obj1)/(Obj1+obj2)'})

% subplot(2,2,3)
% [hb, eb, sc, p, stats] = MS_bar_w_err(this_data(r_idx & ctrl_idx)', this_data(r_idx & (art_idx | bi_idx))', [[ .6 .6 .6];c_ord(1,:)],1,  'ttest2', [1 2]);


% plot the interaction number and events per type as a scatter
subplot(2,2,2)
hold on
scatter(tbl_out.DI_t(r_idx & ctrl_idx), tbl_out.DI_n(r_idx & ctrl_idx), 100, 'MarkerFaceColor', [.6 .6 .6],"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
scatter(tbl_out.DI_t(r_idx & bi_idx), tbl_out.DI_n(r_idx & bi_idx), 100, 'MarkerFaceColor', c_ord(2,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
scatter(tbl_out.DI_t(r_idx & art_idx), tbl_out.DI_n(r_idx & art_idx), 100, 'MarkerFaceColor', c_ord(3,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
scatter(tbl_out.DI_t(r_idx & cheta_idx), tbl_out.DI_n(r_idx & cheta_idx), 100, 'MarkerFaceColor', c_ord(5,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
xline(0); yline(0);
ylabel('DI: n interactions')
xlabel('DI: time')
legend({'Control', 'Bipole', 'ArchT', 'Cheta'}, 'Box','off', 'Location','northwest')
xlim([-1 1]); ylim([-1 1])

subplot(2,2,3)
cla
hold on
yline(0)

% ctrl
idx = find(r_idx & ctrl_idx);
x_vec = sort(MS_randn_range(1,sum(r_idx & ctrl_idx), .8, 1.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', [.6 .6 .6],"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
for ii = 1:length(idx) 
    if contains(tbl_out.Sex(idx(ii)), 'M')
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)))
    end
    disp(tbl_out.Subject(idx(ii)))
end


% bipole
idx = find(r_idx & bi_idx & ~sst_idx);
x_vec = sort(MS_randn_range(1,sum(r_idx & bi_idx & ~sst_idx), 1.8, 2.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', c_ord(2,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
for ii = 1:length(idx) 
    if contains(tbl_out.Sex(idx(ii)), 'M')
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)))
    end
end

idx = find(r_idx & bi_idx & sst_idx);
x_vec = sort(MS_randn_range(1,sum(r_idx & bi_idx & sst_idx), 1.8, 2.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', c_ord(2,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"diamond")
for ii = 1:length(idx) 
    if contains(tbl_out.Sex(idx(ii)), 'M')
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)))
    end
end

% archt
idx = find(r_idx & art_idx & ~sst_idx);
x_vec = sort(MS_randn_range(1,sum(r_idx & art_idx& ~sst_idx), 2.8, 3.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', c_ord(3,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
for ii = 1:length(idx) 
    if contains(tbl_out.Sex(idx(ii)), 'M')
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)))
    end
end

idx = find(r_idx & art_idx & sst_idx);
x_vec = sort(MS_randn_range(1,sum(r_idx & art_idx& sst_idx), 2.8, 3.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', c_ord(3,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"diamond")
for ii = 1:length(idx) 
    if contains(tbl_out.Sex(idx(ii)), 'M')
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)))
    end
end

idx = find(r_idx & cheta_idx & ~sst_idx);
x_vec = sort(MS_randn_range(1,sum(r_idx & cheta_idx & ~sst_idx),3.8, 4.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', c_ord(5,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"o")
for ii = 1:length(idx)
    if contains(tbl_out.Sex(idx(ii)), 'M')
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)))
    end
end

idx = find(r_idx & cheta_idx & sst_idx);
x_vec = sort(MS_randn_range(1,sum(r_idx & cheta_idx & sst_idx),3.8, 4.2));
scatter(x_vec, this_data(idx), 100, 'MarkerFaceColor', c_ord(5,:),"MarkerEdgeColor",[.6 .6 .6], 'Marker',"diamond")
for ii = 1:length(idx)
    if contains(tbl_out.Sex(idx(ii)), 'Male')
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)), 'color', 'r')
    else
        text(x_vec(ii), this_data(idx(ii)), tbl_out.Subject(idx(ii)))
    end
end

set(gca, 'XTick', 1:4, 'XTickLabel', { 'Ctrl - Rec' 'Bipole - Rec'...
    'ArchT - Rec' 'Cheta - Rec'}, 'XTickLabelRotation', 45)

ylim([-1 1]); 
xlim ([-.5 4.5])
ylabel({'Discrimination index'; '(Obj2 - Obj1)/(Obj1+obj2)'})
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TFC 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the LED on table

TFC_tab = readtable('CIHR_TFC_Frame - Sheet1.csv');
for ii = 1:length(TFC_tab.Subject)
    TFC_tab.S_num(ii) = str2double(TFC_tab.Subject{ii}(2:end)); 
end

% protocol

TFC1.baseline = [0 300];
TFC1.tone = [300 320; 492 512; 684 704 ; 876 896; 1068 1088]; 
TFC1.trace = [320 340; 512 532; 704 724; 896 916; 1088 1108]; 
TFC1.shock = [340 342; 532 534; 724 726; 916 918; 1108 1110]; 
TFC1.ITI = [342 492; 534 684; 726 876; 918 1068; 1110 1260];

TFC2.baseline = [0 300];
TFC2.tone = [300 320; 492 512; 684 704; 876 896; 1068 1088];
TFC2.trace = [322 342; 514 534 ;706 726 ;898 918 ;1090 1110]; 
TFC2.ITI = [342 492; 534 684 ;726 876 ;918 1068 ;1110 1260]; 

TFC3.baseline = [0 300];

figure(1010)
clf
subplot(2,1,1)
hold on

blocks = fieldnames(TFC1); 
c_ord = parula(length(blocks)); 

for ii = 1:length(blocks)
    for jj = 1:size(TFC1.(blocks{ii}),1)
        rectangle('Position',[TFC1.(blocks{ii})(jj,1), ii-1, TFC1.(blocks{ii})(jj,2) - TFC1.(blocks{ii})(jj,1), 1], 'FaceColor',c_ord(ii,:));
    end
end
set(gca, 'yTick', .5:length(blocks), 'yTickLabel', blocks)


subplot(2,1,2)
hold on

blocks = fieldnames(TFC2); 
c_ord = parula(length(blocks)); 

for ii = 1:length(blocks)
    for jj = 1:size(TFC2.(blocks{ii}),1)
        rectangle('Position',[TFC2.(blocks{ii})(jj,1), ii-1, TFC2.(blocks{ii})(jj,2) - TFC2.(blocks{ii})(jj,1), 1], 'FaceColor',c_ord(ii,:));
    end
end
set(gca, 'yTick', .5:length(blocks), 'yTickLabel', blocks)
%% initialize

day_dir = dir(cd); 
rm_idx = true(1, length(day_dir)); 
for ii  = 1:length(day_dir); if contains(day_dir(ii).name, 'TFC'); rm_idx(ii) = false; end; end
day_dir(rm_idx) = []; 
f_list = [];

for iD = 1:length(day_dir)

    d_list = dir([day_dir(iD).folder filesep day_dir(iD).name]);
    s_idx = [];
    rm_idx = true(1, length(d_list));
    for ii  = 1:length(d_list); if strcmp(d_list(ii).name(1), 'M'); rm_idx(ii) = false; end; s_idx(ii) = str2double(d_list(ii).name(2:end)); end
    d_list(rm_idx) = []; s_idx(rm_idx) = [];

    [~, s] = sort(s_idx);
    d_list = d_list(s);


    for ii  = length(d_list):-1:1

        p_dir = dir(fullfile([d_list(ii).folder filesep d_list(ii).name], ['**' filesep '*Processed*'])); % Lists all files and folders recursively

        f_list(end+1).folder = p_dir.folder;
    end

    for iF = 1:length(f_list)
        fprintf('%s\n',f_list(iF).folder)
    end


end
%% loop over sessons
out = [];

for iF = 1:length(f_list)

    parts = strsplit(f_list(iF).folder, filesep); 

    s_idx = find(contains(parts, 'TFC_')); 
    
    info = [];
    info.subject = parts{s_idx+2}; 
    info.sess = parts{s_idx+1}; 

    if strcmp(info.sess, 'TFC_A'); info.sess = 'TFC1'; 
    elseif strcmp(info.sess, 'TFC_B'); info.sess = 'TFC2';
    elseif strcmp(info.sess, 'TFC_C'); info.sess = 'TFC3';
    end
    
    if strcmp(info.sess, 'TFC1')
        proto = TFC1;
    elseif strcmp(info.sess, 'TFC2')
        proto = TFC2;
    elseif strcmp(info.sess, 'TFC3')
        proto = TFC3;
    end
    
    % get the table info for the lED on frame.
    this_tab = find(TFC_tab.S_num == str2double(info.subject(2:end)));
    if isempty(this_tab)
        continue
    end
    
    if (~contains(info.subject, {'M5' 'M17'}) && ~isempty(this_tab)) &&...
            (~(contains(info.subject, {'M8'}) && contains(info.sess, {'TFC1'})) && ~isempty(this_tab)) && ...
             (~(contains(info.subject, {'M3'}) && contains(info.sess, {'TFC1'})) && ~isempty(this_tab)) && ...
             (~(strcmp(info.subject(1:end), 'M1') && contains(info.sess, {'TFC1'})) && ~isempty(this_tab)) && ...
             (~(contains(info.subject, {'M7'}) && contains(info.sess, {'TFC1'})) && ~isempty(this_tab)) %+ ~isnan(TFC_tab.(info.sess)(this_tab))) == 0
        
        out.(info.subject).(info.sess) = MS_DLC_score_freezing_dir(f_list(iF).folder,.1,proto, TFC_tab.(info.sess)(this_tab), [MS_parent_dir(d_list(1).folder) filesep 'figs' filesep info.subject '_' info.sess]);
        
    else
        out.(info.subject).(info.sess) = [];
        out.(info.subject).(info.sess).fvec = [];
        out.(info.subject).(info.sess).f_bin = [];
        out.(info.subject).(info.sess).t_bin = [];
        out.(info.subject).(info.sess).TFC = [];
    end
    
    % hold the 60 binned freezing.
    %     if strcmp(info.sess, 'TFC1')
    %         TFC1_out = [TFC1_out, out.(info.subject).(info.sess).f_bin];
    %     elseif strcmp(info.sess, 'TFC2')
    %         TFC2_out = [TFC2_out, out.(info.subject).(info.sess).f_bin];
    %     elseif strcmp(info.sess, 'TFC3')
    %         TFC3_out = [TFC3_out, out.(info.subject).(info.sess).f_bin];
    %     end
    %
end
%

%% collect the outputs

c_ord = MS_linspecer(4); 

s_list = fieldnames(out);
TFC1_out = []; TFC2_out = []; TFC3_out = [];TFC2_bar = [];

TFC1_bar.baseline = [];
TFC1_bar.tone = [];
TFC1_bar.trace = [];
TFC1_bar.ITI = [];

TFC2_bar.baseline = [];
TFC2_bar.tone = [];
TFC2_bar.trace = [];
TFC2_bar.ITI = [];
TFC3_bar.baseline = []; 
Geno = []; Subs = []; Sess = []; 
TFC1_geno = []; TFC2_geno = []; TFC3_geno = []; 
TFC2_sub = []; 

for iSub = 1:length(s_list)
    
    sess_list = fieldnames(out.(s_list{iSub}));
    
    for iSess = 1:length(sess_list)
        if isempty(out.(s_list{iSub}).(sess_list{iSess}).TFC)
            continue
        else
                        this_tab = find(TFC_tab.S_num == str2double(s_list{iSub}(2:end)));

            % hold the 60 binned freezing.
            if strcmp(sess_list{iSess}, 'TFC1')
                TFC1_out = [TFC1_out; out.(s_list{iSub}).(sess_list{iSess}).f_bin];
                t_vec = out.(s_list{iSub}).(sess_list{iSess}).t_bin;
                
                TFC1_bar.baseline(end+1) = mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.baseline);
                TFC1_bar.tone(end+1) =  mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.tone);
                TFC1_bar.trace(end+1) =  mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.trace);
                TFC1_bar.ITI(end+1)=  mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.ITI);

                TFC1_geno(end+1) = TFC_tab.Pox(this_tab); 
                disp([s_list{iSub} ' ' sess_list{iSess}])

            elseif strcmp(sess_list{iSess}, 'TFC2')
                TFC2_out = [TFC2_out; out.(s_list{iSub}).(sess_list{iSess}).f_bin];
                t_vec2 = out.(s_list{iSub}).(sess_list{iSess}).t_bin;
                
                TFC2_bar.baseline(end+1) = mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.baseline);
                TFC2_bar.tone(end+1) =  mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.tone);
                TFC2_bar.trace(end+1) =  mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.trace);
                TFC2_bar.ITI(end+1)=  mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.ITI);
                TFC2_geno(end+1) = TFC_tab.Pox(this_tab);
                TFC2_sub{end+1} = s_list{iSub}; 

            elseif strcmp(sess_list{iSess}, 'TFC3')
                TFC3_out = [TFC3_out; out.(s_list{iSub}).(sess_list{iSess}).f_bin];
                t_vec3 = out.(s_list{iSub}).(sess_list{iSess}).t_bin;
                TFC3_bar.baseline(end+1) = mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.baseline);
                TFC3_geno(end+1) = TFC_tab.Pox(this_tab); 
            end
            Subs{end+1} = s_list{iSub}; 
            Sess{end+1} = sess_list{iSess}; 

            Geno(end+1) = TFC_tab.Pox(this_tab); 
        end
    end
    
end

%%

% simple plots

figure(202)
clf
subplot(3,3,1:2)
hold on
for ii = 1:length(TFC1.tone)
   rectangle('position',  [TFC1.tone(ii,1)/60 0, (TFC1.tone(ii,2) - TFC1.tone(ii,1))/60, 1], 'FaceColor', [.7 .7 .7 .2], 'EdgeColor', 'none', 'FaceAlpha',.4)
   rectangle('position',  [TFC1.trace(ii,1)/60 0, (TFC1.trace(ii,2) - TFC1.trace(ii,1))/60, 1], 'FaceColor', c_ord(3,:), 'EdgeColor', 'none', 'FaceAlpha',.4)
   rectangle('position',  [TFC1.shock(ii,1)/60 0, (TFC1.shock(ii,2) - TFC1.shock(ii,1))/60, 1], 'FaceColor',[([210 240 10]/255) .8], 'EdgeColor', 'none', 'FaceAlpha',.4)

end
p(1) = plot(t_vec/60, mean(TFC1_out(TFC1_geno == 0,:)), 'color', c_ord(1,:), 'linewidth', 3);
errorbar(t_vec/60, mean(TFC1_out(TFC1_geno == 0,:)),MS_SEM_vec(TFC1_out(TFC1_geno == 0,:)), 'color', c_ord(1,:)); 
errorbar(t_vec/60, mean(TFC1_out(TFC1_geno == 0,:)),MS_SEM_vec(TFC1_out(TFC1_geno == 0,:)), 'color', c_ord(1,:)); 

p(2) = plot(t_vec/60, mean(TFC1_out(TFC1_geno == 2,:)), 'color', c_ord(3,:), 'linewidth', 3);
errorbar(t_vec/60, mean(TFC1_out(TFC1_geno == 2,:)),MS_SEM_vec(TFC1_out(TFC1_geno == 2,:)), 'color', c_ord(1,:)); 
errorbar(t_vec/60, mean(TFC1_out(TFC1_geno == 2,:)),MS_SEM_vec(TFC1_out(TFC1_geno == 2,:)), 'color', c_ord(1,:)); 

p(3) = plot(t_vec/60, mean(TFC1_out(TFC1_geno == 1,:)), 'color', c_ord(2,:), 'linewidth', 3);
errorbar(t_vec/60, mean(TFC1_out(TFC1_geno == 1,:)),MS_SEM_vec(TFC1_out(TFC1_geno == 1,:)), 'color', c_ord(2,:)); 
errorbar(t_vec/60, mean(TFC1_out(TFC1_geno == 1,:)),MS_SEM_vec(TFC1_out(TFC1_geno == 1,:)), 'color', c_ord(2,:)); 

p(4) = plot(t_vec/60, mean(TFC1_out(TFC1_geno == 3,:)), 'color', c_ord(4,:), 'linewidth', 3);
errorbar(t_vec/60, mean(TFC1_out(TFC1_geno == 3,:)),MS_SEM_vec(TFC1_out(TFC1_geno == 3,:)), 'color', c_ord(4,:)); 
errorbar(t_vec/60, mean(TFC1_out(TFC1_geno == 3,:)),MS_SEM_vec(TFC1_out(TFC1_geno == 3,:)), 'color', c_ord(4,:)); 


title('Context A: Aquisition'); 
ylabel('%freezing')
xlabel('time (s)')
ylim([0 1])
xlim([0 t_vec(end)/60])

legend([p(1) p(2) p(3) p(4)], {'Ctrl', 'Bipole', 'ArchT', 'Cheta'}, 'Location','southeast', 'Box','off')

subplot(3,3,4:5)
cla
hold on
for ii = 1:length(TFC2.tone)
   rectangle('position',  [TFC2.tone(ii,1)/60 0, (TFC2.tone(ii,2) - TFC2.tone(ii,1))/60, 1], 'FaceColor', [.7 .7 .7 .2], 'EdgeColor', 'none', 'FaceAlpha',.4)
   rectangle('position',  [TFC2.trace(ii,1)/60 0, (TFC2.trace(ii,2) - TFC2.trace(ii,1))/60, 1], 'FaceColor', c_ord(3,:), 'EdgeColor', 'none', 'FaceAlpha',.4)

end

plot(t_vec2/60, mean(TFC2_out(TFC1_geno == 0,:)), 'color', c_ord(1,:), 'linewidth', 3)
errorbar(t_vec2/60, mean(TFC2_out(TFC1_geno == 0,:)),MS_SEM_vec(TFC2_out(TFC2_geno == 0,:)), 'color', c_ord(1,:)); 
errorbar(t_vec2/60, mean(TFC2_out(TFC1_geno == 0,:)),MS_SEM_vec(TFC2_out(TFC2_geno == 0,:)), 'color', c_ord(1,:)); 

plot(t_vec2/60, mean(TFC2_out(TFC2_geno == 1,:)), 'color', c_ord(2,:), 'linewidth', 3)
errorbar(t_vec2/60, mean(TFC2_out(TFC2_geno == 1,:)),MS_SEM_vec(TFC2_out(TFC2_geno == 1,:)), 'color', c_ord(2,:)); 
errorbar(t_vec2/60, mean(TFC2_out(TFC2_geno == 1,:)),MS_SEM_vec(TFC2_out(TFC2_geno == 1,:)), 'color', c_ord(2,:)); 

plot(t_vec2/60, mean(TFC2_out(TFC2_geno == 2,:)), 'color', c_ord(3,:), 'linewidth', 3)
errorbar(t_vec2/60, mean(TFC2_out(TFC2_geno == 2,:)),MS_SEM_vec(TFC2_out(TFC2_geno == 2,:)), 'color', c_ord(3,:)); 
errorbar(t_vec2/60, mean(TFC2_out(TFC2_geno == 2,:)),MS_SEM_vec(TFC2_out(TFC2_geno == 2,:)), 'color', c_ord(3,:)); 


plot(t_vec2/60, mean(TFC2_out(TFC2_geno == 3,:)), 'color', c_ord(4,:), 'linewidth', 3)
errorbar(t_vec2/60, mean(TFC2_out(TFC2_geno == 3,:)),MS_SEM_vec(TFC2_out(TFC2_geno == 3,:)), 'color', c_ord(4,:)); 
errorbar(t_vec2/60, mean(TFC2_out(TFC2_geno == 3,:)),MS_SEM_vec(TFC2_out(TFC2_geno == 3,:)), 'color', c_ord(4,:)); 

title('Context B: Trace Memory'); 
ylabel('%freezing')
xlabel('time (s)')
ylim([0 1])
xlim([0 t_vec2(end)/60])

subplot(3,3,7:8)
hold on
plot(t_vec3/60, mean(TFC3_out(TFC3_geno == 0,:)), 'color', c_ord(1,:), 'linewidth', 3)
errorbar(t_vec3/60, mean(TFC3_out(TFC3_geno == 0,:)),MS_SEM_vec(TFC3_out(TFC2_geno == 0,:)), 'color', c_ord(1,:)); 
errorbar(t_vec3/60, mean(TFC3_out(TFC3_geno == 0,:)),MS_SEM_vec(TFC3_out(TFC2_geno == 0,:)), 'color', c_ord(1,:)); 

plot(t_vec3/60, mean(TFC3_out(TFC3_geno == 1,:)), 'color', c_ord(2,:), 'linewidth', 3)
errorbar(t_vec3/60, mean(TFC3_out(TFC3_geno == 1,:)),MS_SEM_vec(TFC3_out(TFC3_geno == 1,:)), 'color', c_ord(2,:)); 
errorbar(t_vec3/60, mean(TFC3_out(TFC3_geno == 1,:)),MS_SEM_vec(TFC3_out(TFC3_geno == 1,:)), 'color', c_ord(2,:)); 

plot(t_vec3/60, mean(TFC3_out(TFC3_geno == 2,:)), 'color', c_ord(3,:), 'linewidth', 3)
errorbar(t_vec3/60, mean(TFC3_out(TFC3_geno == 2,:)),MS_SEM_vec(TFC3_out(TFC2_geno == 2,:)), 'color', c_ord(3,:)); 
errorbar(t_vec3/60, mean(TFC3_out(TFC3_geno == 2,:)),MS_SEM_vec(TFC3_out(TFC2_geno == 2,:)), 'color', c_ord(3,:)); 

plot(t_vec3/60, mean(TFC3_out(TFC3_geno == 3,:)), 'color', c_ord(4,:), 'linewidth', 3)
errorbar(t_vec3/60, mean(TFC3_out(TFC3_geno == 3,:)),MS_SEM_vec(TFC3_out(TFC3_geno == 3,:)), 'color', c_ord(4,:)); 
errorbar(t_vec3/60, mean(TFC3_out(TFC3_geno == 3,:)),MS_SEM_vec(TFC3_out(TFC3_geno == 3,:)), 'color', c_ord(4,:)); 

xlim([0 t_vec3(end)/60])
title('Context A'': Context Memory'); 
ylabel('%freezing')
xlabel('time (s)')
ylim([0 1])

% bar plots
subplot(3,3,3)
cla
hold on

[hb, eb, sc, p, stats] = MS_bar_w_err(TFC1_bar.baseline(TFC1_geno == 0), TFC1_bar.baseline(TFC1_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [1 2]); 
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
[hb, eb, sc, p, stats] = MS_bar_w_err(TFC1_bar.baseline(TFC1_geno == 2), TFC1_bar.baseline(TFC1_geno == 3), [c_ord(3,:);c_ord(4,:)],1,  'ttest2', [3 4]);
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

[hb, eb, sc, p, stats] = MS_bar_w_err(TFC1_bar.tone(TFC1_geno == 0), TFC1_bar.tone(TFC1_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [6 7]); 
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
[hb, eb, sc, p, stats] = MS_bar_w_err(TFC1_bar.tone(TFC1_geno == 2), TFC1_bar.tone(TFC1_geno == 3), [c_ord(3,:);c_ord(4,:)],1,  'ttest2', [8 9]);
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

[hb, eb, sc, p, stats] = MS_bar_w_err(TFC1_bar.trace(TFC1_geno == 0), TFC1_bar.trace(TFC1_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [11 12]); 
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
[hb, eb, sc, p, stats] = MS_bar_w_err(TFC1_bar.trace(TFC1_geno == 2), TFC1_bar.trace(TFC1_geno == 3), [c_ord(3,:);c_ord(4,:)],1,  'ttest2', [13 14]);
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

[hb, eb, sc, p, stats] = MS_bar_w_err(TFC1_bar.ITI(TFC1_geno == 0), TFC1_bar.ITI(TFC1_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [16 17]); 
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
[hb, eb, sc, p, stats] = MS_bar_w_err(TFC1_bar.ITI(TFC1_geno == 2), TFC1_bar.ITI(TFC1_geno == 3), [c_ord(3,:);c_ord(4,:)],1,  'ttest2', [18 19]);
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';


% scatter(1+sort(MS_randn_range(length(TFC1_bar.tone(TFC1_geno == 0)), 1, -.1, .1)), TFC1_bar.tone(TFC1_geno == 0),25,  c_ord(1,:), 'filled')
% scatter(2+sort(MS_randn_range(length(TFC1_bar.tone(TFC1_geno == 1)), 1, -.1, .1)), TFC1_bar.tone(TFC1_geno == 1),25,  c_ord(2,:), 'filled')
% scatter(3+sort(MS_randn_range(length(TFC1_bar.tone(TFC1_geno == 2)), 1, -.1, .1)), TFC1_bar.tone(TFC1_geno == 2),25,  c_ord(3,:), 'filled')
% scatter(4+sort(MS_randn_range(length(TFC1_bar.tone(TFC1_geno == 3)), 1, -.1, .1)), TFC1_bar.tone(TFC1_geno == 3),25,  c_ord(4,:), 'filled')
% 
% scatter(7+sort(MS_randn_range(length(TFC1_bar.trace(TFC1_geno == 0)), 1, -.1, .1)), TFC1_bar.trace(TFC1_geno == 0),25,  c_ord(1,:), 'filled')
% scatter(8+sort(MS_randn_range(length(TFC1_bar.trace(TFC1_geno == 1)), 1, -.1, .1)), TFC1_bar.trace(TFC1_geno == 1),25,  c_ord(2,:), 'filled')
% scatter(9+sort(MS_randn_range(length(TFC1_bar.trace(TFC1_geno == 2)), 1, -.1, .1)), TFC1_bar.trace(TFC1_geno == 2),25,  c_ord(3,:), 'filled')
% scatter(10+sort(MS_randn_range(length(TFC1_bar.trace(TFC1_geno == 3)), 1, -.1, .1)), TFC1_bar.trace(TFC1_geno == 3),25,  c_ord(4,:), 'filled')
% 
% scatter(13+sort(MS_randn_range(length(TFC1_bar.ITI(TFC1_geno == 0)), 1, -.1, .1)), TFC1_bar.ITI(TFC1_geno == 0),25,  c_ord(1,:), 'filled')
% scatter(14+sort(MS_randn_range(length(TFC1_bar.ITI(TFC1_geno == 1)), 1, -.1, .1)), TFC1_bar.ITI(TFC1_geno == 1),25,  c_ord(2,:), 'filled')
% scatter(15+sort(MS_randn_range(length(TFC1_bar.ITI(TFC1_geno == 2)), 1, -.1, .1)), TFC1_bar.ITI(TFC1_geno == 2),25,  c_ord(3,:), 'filled')
% scatter(16+sort(MS_randn_range(length(TFC1_bar.ITI(TFC1_geno == 3)), 1, -.1, .1)), TFC1_bar.ITI(TFC1_geno == 3),25,  c_ord(4,:), 'filled')
y_l = ylim; 
ylim([y_l(1) 1.1])
set(gca, 'xtick',[2.5 7.5 12.5 17.5], 'XTickLabel', {'baseline', 'tone', 'trace', 'ITI'})
ylabel('Freezing (%)')

subplot(3,3,6)
cla
hold on
[hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.baseline(TFC2_geno == 0), TFC2_bar.baseline(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [1 2]); 
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
[hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.baseline(TFC2_geno == 2), TFC2_bar.baseline(TFC2_geno == 3), [c_ord(3,:);c_ord(4,:)],1,  'ttest2', [3 4]);
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

[hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.tone(TFC2_geno == 0), TFC2_bar.tone(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [6 7]); 
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
[hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.tone(TFC2_geno == 2), TFC2_bar.tone(TFC2_geno == 3), [c_ord(3,:);c_ord(4,:)],1,  'ttest2', [8 9]);
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

[hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.trace(TFC2_geno == 0), TFC2_bar.trace(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [11 12]); 
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
[hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.trace(TFC2_geno == 2), TFC2_bar.trace(TFC2_geno == 3), [c_ord(3,:);c_ord(4,:)],1,  'ttest2', [13 14]);
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';

[hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.ITI(TFC2_geno == 0), TFC2_bar.ITI(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [16 17]); 
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';
% 
[hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.ITI(TFC2_geno == 2), TFC2_bar.ITI(TFC2_geno == 3), [c_ord(3,:);c_ord(4,:)],1,  'ttest2', [18 19]);
hb(1).FaceColor = 'none'; hb(1).EdgeColor = 'k';


% scatter(1+sort(MS_randn_range(length(TFC2_bar.tone(TFC2_geno == 0)), 1, -.1, .1)), TFC2_bar.tone(TFC2_geno == 0),25,  c_ord(1,:), 'filled')
% scatter(2+sort(MS_randn_range(length(TFC2_bar.tone(TFC2_geno == 1)), 1, -.1, .1)), TFC2_bar.tone(TFC2_geno == 1),25,  c_ord(2,:), 'filled')
% scatter(3+sort(MS_randn_range(length(TFC2_bar.tone(TFC2_geno == 2)), 1, -.1, .1)), TFC2_bar.tone(TFC2_geno == 2),25,  c_ord(3,:), 'filled')
% scatter(4+sort(MS_randn_range(length(TFC2_bar.tone(TFC2_geno == 3)), 1, -.1, .1)), TFC2_bar.tone(TFC2_geno == 3),25,  c_ord(4,:), 'filled')
% 
% scatter(1+sort(MS_randn_range(length(TFC2_bar.trace(TFC2_geno == 0)), 1, -.1, .1)), TFC2_bar.trace(TFC2_geno == 0),25,  c_ord(1,:), 'filled')
% scatter(2+sort(MS_randn_range(length(TFC2_bar.trace(TFC2_geno == 1)), 1, -.1, .1)), TFC2_bar.trace(TFC2_geno == 1),25,  c_ord(2,:), 'filled')
% scatter(3+sort(MS_randn_range(length(TFC2_bar.trace(TFC2_geno == 2)), 1, -.1, .1)), TFC2_bar.trace(TFC2_geno == 2),25,  c_ord(3,:), 'filled')
% scatter(4+sort(MS_randn_range(length(TFC2_bar.trace(TFC2_geno == 3)), 1, -.1, .1)), TFC2_bar.trace(TFC2_geno == 3),25,  c_ord(4,:), 'filled')
% 
% scatter(1+sort(MS_randn_range(length(TFC2_bar.ITI(TFC2_geno == 0)), 1, -.1, .1)), TFC2_bar.ITI(TFC2_geno == 0),25,  c_ord(1,:), 'filled')
% scatter(2+sort(MS_randn_range(length(TFC2_bar.ITI(TFC2_geno == 1)), 1, -.1, .1)), TFC2_bar.ITI(TFC2_geno == 1),25,  c_ord(2,:), 'filled')
% scatter(3+sort(MS_randn_range(length(TFC2_bar.ITI(TFC2_geno == 2)), 1, -.1, .1)), TFC2_bar.ITI(TFC2_geno == 2),25,  c_ord(3,:), 'filled')
% scatter(4+sort(MS_randn_range(length(TFC2_bar.ITI(TFC2_geno == 3)), 1, -.1, .1)), TFC2_bar.ITI(TFC2_geno == 3),25,  c_ord(4,:), 'filled')
y_l = ylim; 
ylim([y_l(1) 1.1])
set(gca, 'xtick',[2.5 7.5 12.5 17.5], 'XTickLabel', {'baseline', 'tone', 'trace', 'ITI'})
ylabel('Freezing (%)')

% 
subplot(3,3,9)
cla
hold on
[hb, eb, sc, p, stats] = MS_bar_w_err(TFC3_bar.baseline(TFC3_geno == 0), TFC3_bar.baseline(TFC3_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [1 2]);
% 
hb(1).FaceColor = 'none';
hb(1).EdgeColor = 'k';
% 
[hb, eb, sc, p, stats] = MS_bar_w_err(TFC3_bar.baseline(TFC3_geno == 2), TFC3_bar.baseline(TFC3_geno == 3), [c_ord(3,:);c_ord(4,:)],1,  'ttest2', [3 4]);
hb(1).FaceColor = 'none';
hb(1).EdgeColor = 'k';
% 
% [hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.trace(TFC2_geno == 0), TFC2_bar.trace(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [7 8]);
% hb(1).FaceColor = 'none';
% hb(1).EdgeColor = 'k';
% 
% [hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.ITI(TFC2_geno == 0), TFC2_bar.ITI(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [10 11]);
% hb(1).FaceColor = 'none';
% hb(1).EdgeColor = 'k';
% 
% 
% [hb, eb, sc, p, stats] = MS_bar_w_err(TFC3_bar.baseline(TFC3_geno == 0), TFC3_bar.baseline(TFC3_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [14 15]);
% hb(1).FaceColor = 'none';
% hb(1).EdgeColor = 'k';

y_l = ylim; 
ylim([y_l(1) 1.1])
set(gca, 'xtick',[1:4], 'XTickLabel', {'Control', 'Bipole', 'ArchT', 'Cheta'})
ylabel('Freezing (%)')

% text(1.5, y_l(2)*1.1, 'Tone', 'HorizontalAlignment','center','VerticalAlignment','top')

% text(5.5,y_l(2)*1.1, 'Trace', 'HorizontalAlignment','center','VerticalAlignment','top')

% set(gca, 'xtick', [1:2 5:6], 'xticklabels', {'Tau -', 'Tau +', 'Tau -', 'Tau +'}, 'XTickLabelRotation', 45)
% set(gca, 'xtick', [1.5 4.5 7.5 10.5 14.5], 'XTickLabel', {'Base', 'Tone', 'Trace', 'ITI', 'Context'})

% y_l = ylim; 



% title('Context B')
% ylabel('Trace freezing (%)')

% 
% subplot(2,3,5:6)
% hold on
% [hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.baseline(TFC2_geno == 0), TFC2_bar.baseline(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [1 2]);
% 
% hb(1).FaceColor = 'none';
% hb(1).EdgeColor = 'k';
% 
% [hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.tone(TFC2_geno == 0), TFC2_bar.tone(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [4 5]);
% hb(1).FaceColor = 'none';
% hb(1).EdgeColor = 'k';
% 
% [hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.trace(TFC2_geno == 0), TFC2_bar.trace(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [7 8]);
% hb(1).FaceColor = 'none';
% hb(1).EdgeColor = 'k';
% 
% [hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.ITI(TFC2_geno == 0), TFC2_bar.ITI(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [10 11]);
% hb(1).FaceColor = 'none';
% hb(1).EdgeColor = 'k';
% 
% 
% [hb, eb, sc, p, stats] = MS_bar_w_err(TFC3_bar.baseline(TFC3_geno == 0), TFC3_bar.baseline(TFC3_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [14 15]);
% hb(1).FaceColor = 'none';
% hb(1).EdgeColor = 'k';
% 
% y_l = ylim; 
% ylim([y_l(1) 1.1])
% 
% % text(1.5, y_l(2)*1.1, 'Tone', 'HorizontalAlignment','center','VerticalAlignment','top')
% 
% % text(5.5,y_l(2)*1.1, 'Trace', 'HorizontalAlignment','center','VerticalAlignment','top')
% 
% % set(gca, 'xtick', [1:2 5:6], 'xticklabels', {'Tau -', 'Tau +', 'Tau -', 'Tau +'}, 'XTickLabelRotation', 45)
% set(gca, 'xtick', [1.5 4.5 7.5 10.5 14.5], 'XTickLabel', {'Base', 'Tone', 'Trace', 'ITI', 'Context'})


% % scatter(1+sort(MS_randn_range(length(TFC3_bar.baseline(TFC3_geno == 0)), 1, -.1, .1)), TFC3_bar.baseline(TFC3_geno == 0),25,  c_ord(1,:), 'filled')
% % scatter(2+sort(MS_randn_range(length(TFC3_bar.baseline(TFC3_geno == 1)), 1, -.1, .1)), TFC3_bar.baseline(TFC3_geno == 1),25,  c_ord(2,:), 'filled')
% 
% [hb, eb, sc, p, stats] = MS_bar_w_err(TFC3_bar.baseline(TFC3_geno == 0), TFC3_bar.baseline(TFC3_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2');
% 
% set(gca, 'xtick', 1:2, 'xticklabels', {'Tau -', 'Tau +'}, 'XTickLabelRotation', 45)
% ylim([y_l(1) y_l(2)*1.1])
% hb(1).FaceColor = 'none';
% hb(1).EdgeColor = 'k';
% xlim([0 12])
% 
% ylim(y_l)
% title('Context A Re-exposure')
% ylabel('Freezing (%)')

