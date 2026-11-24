%% Master_Fig1_hypno_summary


% collects all of the hypnograms and computes basic stats and summary
% plots;

hypno_dir = '/Users/ecar/Williams Lab Dropbox/Eric Carmichael/JisooProject2020/2020_Results_aftercutting/Across_episodes/Inter';

inter_dir = '/Users/ecar/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter/Hypno'; 

fig_dir = 'C:\Users\ecar\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\PC9\checks';

if ~exist(inter_dir, "dir"); mkdir(inter_dir); end
%% copy all of the hypnograms (only need to run once)

filelist = dir(fullfile(hypno_dir, ['**' filesep '*Hypno.mat*']));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list

% remove pREM
rm_idx = []; 
for ii = length(filelist):-1:1
    if contains(filelist(ii).folder, 'pREM')
        rm_idx(ii) = 1;
    else
        rm_idx(ii) = 0;
    end
end

filelist(logical(rm_idx)) = []; 

% copy the hypno but named after the session. 

for ii  = 1:length(filelist)
    pv_idx = strfind(lower(filelist(ii).folder), 'pv'); 
    s_name = filelist(ii).folder(pv_idx(end):end); 

    copyfile([filelist(ii).folder filesep filelist(ii).name], [inter_dir filesep 'Hypno_' upper(s_name) '.mat'])

    
end
%%%%%%%%%%%%%%%%
%% Load each hypno and grab all the data


cd(inter_dir)
h_list = dir('Hypno*'); 
all_hypno_ratio = []; 
all_hypno_time = [];

for ii = length(h_list):-1:1
    load(h_list(ii).name)
    Hypno = out.hypno; 
    t_dur = Hypno.tvec(end) - Hypno.tvec(1); 
    w = sum(Hypno.data == 1);
    sws = sum(Hypno.data == 2);
    rem = sum(Hypno.data == 3);
    q = sum(Hypno.data == 4);
    t = sum(Hypno.data == 5); 

    all_hypno_ratio(ii,:) = [w, sws, rem, q, t]./length(Hypno.data); 
    all_hypno_time(ii,:)  = all_hypno_ratio(ii,:).*(t_dur/60); 
    parts = strsplit(h_list(ii).name, '_'); 
    all_hypno_sess{ii} = parts{end}(1:end-4); 
clear Hypno

end

%% split into session types
s_type = unique(all_hypno_sess); 

for ii = length(s_type):-1:1
    this_idx = contains(all_hypno_sess, s_type{ii}); 

    sess_rem_time{ii} = all_hypno_time(this_idx,3)'; 
    sess_rem_ratio{ii} = all_hypno_ratio(this_idx,3)'; 
    sess_sws_time{ii} = all_hypno_time(this_idx,2)'; 
    sess_sws_ratio{ii} = all_hypno_ratio(this_idx,2)'; 

    sess_wake_time{ii} = all_hypno_time(this_idx,2)'; 
    sess_wake_ratio{ii} = all_hypno_ratio(this_idx,2)'; 

end



%%
f_pos = [1 6 6.25 3.4]; 
figure(102)
clf
set(gcf,'Units','inch','OuterPosition',f_pos);


subplot(2,4,1); cla
fprintf('<strong>Hypno Ratio</strong>: \n')
[h, eb, sc, p, stats] = MS_bar_w_err3(all_hypno_ratio(:,1)*100,all_hypno_ratio(:,3)*100,all_hypno_ratio(:,2)*100,[ hex2rgb('#808080'); hex2rgb('#437F97'); hex2rgb('#849323')] , 1, 'ranova', 1:3);
eb.LineWidth = .5; %eb.Color = 'k'; eb.LineStyle = "--"; 
h.LineWidth = .8; h.EdgeColor = "none";
sc{1}.SizeData = 5; sc{2}.SizeData = 5; sc{3}.SizeData = 5; 
sc{1}.MarkerFaceColor = hex2rgb('#808080'); sc{2}.MarkerFaceColor = hex2rgb('#437F97'); sc{3}.MarkerFaceColor = hex2rgb('#849323'); 
sc{1}.MarkerEdgeColor = hex2rgb('#333333'); sc{2}.MarkerEdgeColor = hex2rgb('#333333'); sc{3}.MarkerEdgeColor = hex2rgb('#333333'); 

xlim([.5 3.5])
ylabel('Mean % of time in sleep state')
set(gca, 'xticklabel', {'wake', 'REM', 'SWS'}, 'XTickLabelRotation', 0, 'fontsize', 7);

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength',get(gca, 'TickLength')*2)


% time
subplot(2,4,5); cla
fprintf('<strong>Hypno Ratio</strong>: \n')
[h, eb, sc, p, stats] = MS_bar_w_err3(all_hypno_time(:,1),all_hypno_time(:,3),all_hypno_time(:,2),[ hex2rgb('#808080'); hex2rgb('#437F97'); hex2rgb('#849323')] , 1, 'ranova', 1:3);
eb.LineWidth = .5; %eb.Color = 'k'; eb.LineStyle = "--"; 
h.LineWidth = .8; h.EdgeColor = "none";
sc{1}.SizeData = 5; sc{2}.SizeData = 5; sc{3}.SizeData = 5; 
sc{1}.MarkerFaceColor = hex2rgb('#808080'); sc{2}.MarkerFaceColor = hex2rgb('#437F97'); sc{3}.MarkerFaceColor = hex2rgb('#849323'); 
sc{1}.MarkerEdgeColor = hex2rgb('#333333'); sc{2}.MarkerEdgeColor = hex2rgb('#333333'); sc{3}.MarkerEdgeColor = hex2rgb('#333333'); 

xlim([.5 3.5])
ylabel('Mean time in sleep state')
set(gca, 'xticklabel', {'wake', 'REM', 'SWS'}, 'XTickLabelRotation', 0, 'fontsize', 7);

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength',get(gca, 'TickLength')*2)


%% get the times per session time

% SWS
grp_data = [sess_sws_ratio{4}*100,sess_sws_ratio{5}*100,sess_sws_ratio{1}*100]; 
grp_idx = [ones(size(sess_sws_ratio{4})),ones(size(sess_sws_ratio{5}))+1, ones(size(sess_sws_ratio{1}))+2]; 
grps = unique(grp_idx); 

d_mat = NaN(length(grps), max([sum(grp_idx==1), sum(grp_idx==2), sum(grp_idx==3)]));

    d_mat(1,1:sum(grp_idx==1)) = grp_data(grp_idx == grps(1));
    d_mat(2,1:sum(grp_idx==2)) = grp_data(grp_idx == grps(2));
    d_mat(3,1:sum(grp_idx==3)) = grp_data(grp_idx == grps(3));


subplot(2,4,2); cla
fprintf('<strong>Hypno Ratio</strong>: \n')
[h, eb, sc, p, stats] = MS_bar_w_err3(d_mat(1,:),d_mat(2,:), d_mat(3,:),[ hex2rgb('#808080'); hex2rgb('#437F97'); hex2rgb('#849323')] , 1, 'anova1', 1:3);
eb.LineWidth = .5; %eb.Color = 'k'; eb.LineStyle = "--"; 
h.LineWidth = .8; h.EdgeColor = "none";
sc{1}.SizeData = 5; sc{2}.SizeData = 5; sc{3}.SizeData = 5; 
sc{1}.MarkerFaceColor = hex2rgb('#808080'); sc{2}.MarkerFaceColor = hex2rgb('#437F97'); sc{3}.MarkerFaceColor = hex2rgb('#849323'); 
sc{1}.MarkerEdgeColor = hex2rgb('#333333'); sc{2}.MarkerEdgeColor = hex2rgb('#333333'); sc{3}.MarkerEdgeColor = hex2rgb('#333333'); 

xlim([.5 3.5])
ylabel('Mean %SWS')
set(gca, 'xticklabel', {'Nov', 'Fam', 'Anx'}, 'XTickLabelRotation', 0, 'fontsize', 7);

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength',get(gca, 'TickLength')*2)


grp_data = [sess_sws_time{4},sess_sws_time{5},sess_sws_time{1}]; 
grp_idx = [ones(size(sess_sws_ratio{4})),ones(size(sess_sws_ratio{5}))+1, ones(size(sess_sws_ratio{1}))+2]; 
grps = unique(grp_idx); 

d_mat = NaN(length(grps), max([sum(grp_idx==1), sum(grp_idx==2), sum(grp_idx==3)]));
    d_mat(1,1:sum(grp_idx==1)) = grp_data(grp_idx == grps(1));
    d_mat(2,1:sum(grp_idx==2)) = grp_data(grp_idx == grps(2));
    d_mat(3,1:sum(grp_idx==3)) = grp_data(grp_idx == grps(3));

subplot(2,4,6); cla
fprintf('<strong>Hypno Ratio</strong>: \n')
[h, eb, sc, p, stats] = MS_bar_w_err3(d_mat(1,:),d_mat(2,:), d_mat(3,:),[ hex2rgb('#808080'); hex2rgb('#437F97'); hex2rgb('#849323')] , 1, 'anova1', 1:3);
eb.LineWidth = .5; %eb.Color = 'k'; eb.LineStyle = "--"; 
h.LineWidth = .8; h.EdgeColor = "none";
sc{1}.SizeData = 5; sc{2}.SizeData = 5; sc{3}.SizeData = 5; 
sc{1}.MarkerFaceColor = hex2rgb('#808080'); sc{2}.MarkerFaceColor = hex2rgb('#437F97'); sc{3}.MarkerFaceColor = hex2rgb('#849323'); 
sc{1}.MarkerEdgeColor = hex2rgb('#333333'); sc{2}.MarkerEdgeColor = hex2rgb('#333333'); sc{3}.MarkerEdgeColor = hex2rgb('#333333'); 

xlim([.5 3.5])
ylabel('Mean time in SWS (min)')
set(gca, 'xticklabel', {'Nov', 'Fam', 'Anx'}, 'XTickLabelRotation', 0, 'fontsize', 7);

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength',get(gca, 'TickLength')*2)




% REM  
grp_data = [sess_rem_ratio{4}*100,sess_rem_ratio{5}*100,sess_rem_ratio{1}*100]; 
grp_idx = [ones(size(sess_rem_ratio{4})),ones(size(sess_rem_ratio{5}))+1, ones(size(sess_rem_ratio{1}))+2]; 
grps = unique(grp_idx); 

d_mat = NaN(length(grps), max([sum(grp_idx==1), sum(grp_idx==2), sum(grp_idx==3)]));

    d_mat(1,1:sum(grp_idx==1)) = grp_data(grp_idx == grps(1));
    d_mat(2,1:sum(grp_idx==2)) = grp_data(grp_idx == grps(2));
    d_mat(3,1:sum(grp_idx==3)) = grp_data(grp_idx == grps(3));


subplot(2,4,3); cla
fprintf('<strong>Hypno Ratio</strong>: \n')
[h, eb, sc, p, stats] = MS_bar_w_err3(d_mat(1,:),d_mat(2,:), d_mat(3,:),[ hex2rgb('#808080'); hex2rgb('#437F97'); hex2rgb('#849323')] , 1, 'anova1', 1:3);
eb.LineWidth = .5; %eb.Color = 'k'; eb.LineStyle = "--"; 
h.LineWidth = .8; h.EdgeColor = "none";
sc{1}.SizeData = 5; sc{2}.SizeData = 5; sc{3}.SizeData = 5; 
sc{1}.MarkerFaceColor = hex2rgb('#808080'); sc{2}.MarkerFaceColor = hex2rgb('#437F97'); sc{3}.MarkerFaceColor = hex2rgb('#849323'); 
sc{1}.MarkerEdgeColor = hex2rgb('#333333'); sc{2}.MarkerEdgeColor = hex2rgb('#333333'); sc{3}.MarkerEdgeColor = hex2rgb('#333333'); 

xlim([.5 3.5])
ylabel('Mean %REM')
set(gca, 'xticklabel', {'Nov', 'Fam', 'Anx'}, 'XTickLabelRotation', 0, 'fontsize', 7);

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength',get(gca, 'TickLength')*2)



grp_data = [sess_rem_time{4},sess_rem_time{5},sess_rem_time{1}]; 
grp_idx = [ones(size(sess_rem_ratio{4})),ones(size(sess_rem_ratio{5}))+1, ones(size(sess_rem_ratio{1}))+2]; 
grps = unique(grp_idx); 

d_mat = NaN(length(grps), max([sum(grp_idx==1), sum(grp_idx==2), sum(grp_idx==3)]));

    d_mat(1,1:sum(grp_idx==1)) = grp_data(grp_idx == grps(1));
    d_mat(2,1:sum(grp_idx==2)) = grp_data(grp_idx == grps(2));
    d_mat(3,1:sum(grp_idx==3)) = grp_data(grp_idx == grps(3));


subplot(2,4,7); cla
fprintf('<strong>Hypno Ratio</strong>: \n')
[h, eb, sc, p, stats] = MS_bar_w_err3(d_mat(1,:),d_mat(2,:), d_mat(3,:),[ hex2rgb('#808080'); hex2rgb('#437F97'); hex2rgb('#849323')] , 1, 'anova1', 1:3);
eb.LineWidth = .5; %eb.Color = 'k'; eb.LineStyle = "--"; 
h.LineWidth = .8; h.EdgeColor = "none";
sc{1}.SizeData = 5; sc{2}.SizeData = 5; sc{3}.SizeData = 5; 
sc{1}.MarkerFaceColor = hex2rgb('#808080'); sc{2}.MarkerFaceColor = hex2rgb('#437F97'); sc{3}.MarkerFaceColor = hex2rgb('#849323'); 
sc{1}.MarkerEdgeColor = hex2rgb('#333333'); sc{2}.MarkerEdgeColor = hex2rgb('#333333'); sc{3}.MarkerEdgeColor = hex2rgb('#333333'); 

xlim([.5 3.5])
ylabel('Mean time in REM (min)')
set(gca, 'xticklabel', {'Nov', 'Fam', 'Anx'}, 'XTickLabelRotation', 0, 'fontsize', 7);

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength',get(gca, 'TickLength')*2)


% wake  
grp_data = [sess_rem_ratio{4}*100,sess_rem_ratio{5}*100,sess_rem_ratio{1}*100]; 
grp_idx = [ones(size(sess_rem_ratio{4})),ones(size(sess_rem_ratio{5}))+1, ones(size(sess_rem_ratio{1}))+2]; 
grps = unique(grp_idx); 

d_mat = NaN(length(grps), max([sum(grp_idx==1), sum(grp_idx==2), sum(grp_idx==3)]));

    d_mat(1,1:sum(grp_idx==1)) = grp_data(grp_idx == grps(1));
    d_mat(2,1:sum(grp_idx==2)) = grp_data(grp_idx == grps(2));
    d_mat(3,1:sum(grp_idx==3)) = grp_data(grp_idx == grps(3));


subplot(2,4,3); cla
fprintf('<strong>Hypno Ratio</strong>: \n')
[h, eb, sc, p, stats] = MS_bar_w_err3(d_mat(1,:),d_mat(2,:), d_mat(3,:),[ hex2rgb('#808080'); hex2rgb('#437F97'); hex2rgb('#849323')] , 1, 'anova1', 1:3);
eb.LineWidth = .5; %eb.Color = 'k'; eb.LineStyle = "--"; 
h.LineWidth = .8; h.EdgeColor = "none";
sc{1}.SizeData = 5; sc{2}.SizeData = 5; sc{3}.SizeData = 5; 
sc{1}.MarkerFaceColor = hex2rgb('#808080'); sc{2}.MarkerFaceColor = hex2rgb('#437F97'); sc{3}.MarkerFaceColor = hex2rgb('#849323'); 
sc{1}.MarkerEdgeColor = hex2rgb('#333333'); sc{2}.MarkerEdgeColor = hex2rgb('#333333'); sc{3}.MarkerEdgeColor = hex2rgb('#333333'); 

xlim([.5 3.5])
ylabel('Mean %REM')
set(gca, 'xticklabel', {'Nov', 'Fam', 'Anx'}, 'XTickLabelRotation', 0, 'fontsize', 7);

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength',get(gca, 'TickLength')*2)



grp_data = [sess_rem_time{4},sess_rem_time{5},sess_rem_time{1}]; 
grp_idx = [ones(size(sess_rem_ratio{4})),ones(size(sess_rem_ratio{5}))+1, ones(size(sess_rem_ratio{1}))+2]; 
grps = unique(grp_idx); 

d_mat = NaN(length(grps), max([sum(grp_idx==1), sum(grp_idx==2), sum(grp_idx==3)]));

    d_mat(1,1:sum(grp_idx==1)) = grp_data(grp_idx == grps(1));
    d_mat(2,1:sum(grp_idx==2)) = grp_data(grp_idx == grps(2));
    d_mat(3,1:sum(grp_idx==3)) = grp_data(grp_idx == grps(3));


subplot(2,4,7); cla
fprintf('<strong>Hypno Ratio</strong>: \n')
[h, eb, sc, p, stats] = MS_bar_w_err3(d_mat(1,:),d_mat(2,:), d_mat(3,:),[ hex2rgb('#808080'); hex2rgb('#437F97'); hex2rgb('#849323')] , 1, 'anova1', 1:3);
eb.LineWidth = .5; %eb.Color = 'k'; eb.LineStyle = "--"; 
h.LineWidth = .8; h.EdgeColor = "none";
sc{1}.SizeData = 5; sc{2}.SizeData = 5; sc{3}.SizeData = 5; 
sc{1}.MarkerFaceColor = hex2rgb('#808080'); sc{2}.MarkerFaceColor = hex2rgb('#437F97'); sc{3}.MarkerFaceColor = hex2rgb('#849323'); 
sc{1}.MarkerEdgeColor = hex2rgb('#333333'); sc{2}.MarkerEdgeColor = hex2rgb('#333333'); sc{3}.MarkerEdgeColor = hex2rgb('#333333'); 

xlim([.5 3.5])
ylabel('Mean time in REM (min)')
set(gca, 'xticklabel', {'Nov', 'Fam', 'Anx'}, 'XTickLabelRotation', 0, 'fontsize', 7);

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength',get(gca, 'TickLength')*2)

%% save
print("-bestfit",[fig_dir filesep 'fig1_hypno'], '-dpdf', "-vector")