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

    t_dur = Hypno.tvec(end) - Hypno.tvec(1); 
    w = sum(Hypno.data == 1);
    sws = sum(Hypno.data == 2);
    rem = sum(Hypno.data == 3);
    q = sum(Hypno.data == 4);
    t = sum(Hypno.data == 5); 

    all_hypno_ratio(ii,:) = [w, sws, rem, q, t]./length(Hypno.data); 
    all_hypno_time(ii,:)  = all_hypno_ratio(ii,:).*(t_dur/60); 
clear Hypno

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
% subplot(2,4,4); cla
% [h, eb, sc, p, stats] = MS_bar_w_err3(all_hypno_time(:,1),all_hypno_ratio(:,2),all_hypno_ratio(:,3),viridis(3), 1, 'anova1', 1:3); 
% sc{1}.SizeData = 10; sc{2}.SizeData = 10; sc{3}.SizeData = 10; 
% ylabel('Mean time in hypno state (min)')
% set(gca, 'xticklabel', {'wake', 'SWS', 'REM'})


%% save
print("-bestfit",[fig_dir filesep 'fig1_hypno'], '-dpdf', "-vector")