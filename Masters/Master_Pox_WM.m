%% Master collector script for Pox WM data




%% load the learning rates for each cohort. 

all = []; 

all.Jan25 = readtable('POX_WM - Jan25_Cohort.csv'); 
all.Jul24 = readtable('POX_WM - Jul24_Cohort.csv'); 






%% loop over cohorts and grab all the subjects

data.day_mean_ctrl = []; 
data.day_mean_tau = []; 
data.cross_ctrl = []; 
data.cross_tau = []; 
data.day_first_ctrl = [];
data.day_first_tau = [];


c_ord = MS_linspecer(4); 

coh = fieldnames(all); 

for iC = 1:length(coh)

    subs = all.(coh{iC}).ID; 
    keep_idx = contains(subs, 'Pox'); 

    this_tab = all.(coh{iC}); 
    this_tab(~keep_idx,:) = []; 

    this_data = table2array(this_tab(:,4:23));
    this_day_mean = []; 
    for ii = size(this_data,1):-1:1
            this_day_mean(ii,1:5) = [mean(this_data(ii,1:4)) mean(this_data(ii,5:8)), mean(this_data(ii,9:12)), mean(this_data(ii,13:16)), mean(this_data(ii,17:20))];
    end
    
    this_day_mean_first = []; 
    for ii = size(this_data,1):-1:1
            this_day_mean_first(ii,1:5) = [mean(this_data(ii,1:2)) mean(this_data(ii,5:6)), mean(this_data(ii,9:10)), mean(this_data(ii,13:15)), mean(this_data(ii,17:18))];
    end

     

    g_temp = table2array(this_tab(:,2)); 
    
    keep_idx = ~isnan(g_temp); 
    
%     g_temp > 0

    data.(coh{iC}).sub = table2cell(this_tab(keep_idx,1));
    
    data.(coh{iC}).geno= logical(g_temp(keep_idx) >0); 

    data.(coh{iC}).data =  this_data(keep_idx,:); 

    data.(coh{iC}).day_mean =  this_day_mean(keep_idx,:); 
    
        data.(coh{iC}).day_first =  this_day_mean_first(keep_idx,:); 

    data.(coh{iC}).day_mean_ctrl =  data.(coh{iC}).day_mean(data.(coh{iC}).geno == 0,:); 

    data.(coh{iC}).day_mean_tau =  data.(coh{iC}).day_mean(data.(coh{iC}).geno == 1,:);
    
    data.(coh{iC}).day_first_ctrl =  data.(coh{iC}).day_first(data.(coh{iC}).geno == 0,:); 

    data.(coh{iC}).day_first_tau =  data.(coh{iC}).day_first(data.(coh{iC}).geno == 1,:);

    data.(coh{iC}).cross = table2array(this_tab(keep_idx,24)); 


% collect data across cohorts

data.day_mean_ctrl = [data.day_mean_ctrl; data.(coh{iC}).day_mean_ctrl]; 
data.day_mean_tau = [data.day_mean_tau; data.(coh{iC}).day_mean_tau ]; 

data.day_first_ctrl = [data.day_first_ctrl; data.(coh{iC}).day_first_ctrl]; 
data.day_first_tau = [data.day_first_tau; data.(coh{iC}).day_first_tau ]; 

data.cross_ctrl = [data.cross_ctrl; data.(coh{iC}).cross(data.(coh{iC}).geno == 0)]; 
data.cross_tau = [data.cross_tau; data.(coh{iC}).cross(data.(coh{iC}).geno == 1)]; 


    figure(iC)
    clf
    title(['Cohort:' num2str(coh{iC})])
    subplot(1, 4, 1:3)
    cla
    hold on

    % for ii = 1:length(data.(coh{iC}).sub)
    %     plot(1:5, this_day_mean(ii,:))
    % 
    % end    
    plot(1:5, mean(data.(coh{iC}).day_mean_ctrl), 'color', c_ord(1,:), 'LineWidth',2)
    errorb(1:5, mean(data.(coh{iC}).day_mean_ctrl),MS_SEM_vec(data.(coh{iC}).day_mean_ctrl), 'color', c_ord(1,:), 'LineWidth',2)
    
    plot(1:5, mean(data.(coh{iC}).day_mean_tau), 'color', c_ord(2,:), 'LineWidth',2)
    errorb(1:5, mean(data.(coh{iC}).day_mean_tau),MS_SEM_vec(data.(coh{iC}).day_mean_tau), 'color', c_ord(2,:), 'LineWidth',2)


%     plot(1:5, mean(this_day_mean), 'k', 'LineWidth',2)
    % leg = data.(coh{iC}).sub; 
    % leg{end+1} = 'mean';
    legend({'Ctrl', 'Tau', 'All'}, 'Box','off')
    ylabel('Mean escape latency')

    set(gca,'xtick', 1:5, 'XTickLabel', {'Day 1' 'Day2 ' 'Day 3' 'Day 4' 'Day 5'}, 'XTickLabelRotation', 45)

    subplot(1,4,4)
        
    b = MS_bar_w_err(data.(coh{iC}).cross(data.(coh{iC}).geno == 0),data.(coh{iC}).cross(data.(coh{iC}).geno == 1),c_ord(1:2,:), 1, 'ttest2' );
    b(1).FaceColor = 'none';
    b(1).EdgeColor = 'k';

    scatter(1+sort(MS_randn_range(length(data.(coh{iC}).cross(data.(coh{iC}).geno == 0)), 1, -.1, .1)), data.(coh{iC}).cross(data.(coh{iC}).geno == 0),25,  c_ord(1,:), 'filled')
    scatter(2+sort(MS_randn_range(length(data.(coh{iC}).cross(data.(coh{iC}).geno == 1)), 1, -.1, .1)), data.(coh{iC}).cross(data.(coh{iC}).geno == 1),25,  c_ord(2,:), 'filled')

    ylabel('Annulus crossings')
    set(gca, 'xtick', 1:2, 'xticklabel', {'Tau -', 'Tau +'}, 'XTickLabelRotation', 45)

end


%% Combined data

   figure(101)
    clf
    title('All Data')
    subplot(1, 5, 1:2)
    cla
    hold on


    plot(1:5, mean(data.day_mean_ctrl), 'color', c_ord(1,:), 'LineWidth',2)
    errorb(1:5, mean(data.day_mean_ctrl),MS_SEM_vec(data.day_mean_ctrl), 'color', c_ord(1,:), 'LineWidth',2)
    
    plot(1:5, mean(data.day_mean_tau), 'color', c_ord(2,:), 'LineWidth',2)
    errorb(1:5, mean(data.day_mean_tau),MS_SEM_vec(data.day_mean_tau), 'color', c_ord(2,:), 'LineWidth',2)


%     plot(1:5, mean(this_day_mean), 'k', 'LineWidth',2)
    % leg = data.(coh{iC}).sub; 
    % leg{end+1} = 'mean';
    legend({'Tau -', 'Tau +', 'All'}, 'Box','off')
   ylabel({'Escape latency (s)'})
    title('All trials')
    ylim([0 40])

    set(gca,'xtick', 1:5, 'XTickLabel', {'Day 1' 'Day2 ' 'Day 3' 'Day 4' 'Day 5'}, 'XTickLabelRotation', 45)
    
    
    % same thing but only first two trials
        subplot(1, 5, 3:4)
    cla
    hold on


    plot(1:5, mean(data.day_first_ctrl), 'color', c_ord(1,:), 'LineWidth',2)
    errorb(1:5, mean(data.day_first_ctrl),MS_SEM_vec(data.day_first_ctrl), 'color', c_ord(1,:), 'LineWidth',2)
    
    plot(1:5, mean(data.day_first_tau), 'color', c_ord(2,:), 'LineWidth',2)
    errorb(1:5, mean(data.day_first_tau),MS_SEM_vec(data.day_first_tau), 'color', c_ord(2,:), 'LineWidth',2)

    
    % anova 
   
%     plot(1:5, mean(this_day_mean), 'k', 'LineWidth',2)
    % leg = data.(coh{iC}).sub; 
    % leg{end+1} = 'mean';
    legend({'Tau -', 'Tau +', 'All'}, 'Box','off')
    ylabel({'Escape latency (s)'})
    title('First two trials')
    ylim([0 40])

    set(gca,'xtick', 1:5, 'XTickLabel', {'Day 1' 'Day2 ' 'Day 3' 'Day 4' 'Day 5'}, 'XTickLabelRotation', 45)


    subplot(1,5,5)
        
    [b, ~, ~, p, stats] = MS_bar_w_err(data.cross_ctrl,data.cross_tau,c_ord(1:2,:), 1, 'ttest2');
    b(1).FaceColor = 'none';
    b(1).EdgeColor = 'k';
    fprintf('Annulus crossing ttest: df: %2.0f  tstat: %2.2f p: %2.2f\n', stats.df, stats.tstat, p)

    scatter(1+sort(MS_randn_range(length(data.cross_ctrl), 1, -.1, .1)), data.cross_ctrl,25,  c_ord(1,:), 'filled')
    scatter(2+sort(MS_randn_range(length(data.cross_tau), 1, -.1, .1)), data.cross_tau,25,  c_ord(2,:), 'filled')

    ylabel('Annulus crossings')
    set(gca, 'xtick', 1:2, 'xticklabel', {'Tau -', 'Tau +'}, 'XTickLabelRotation', 45)
    
    
    cfg_fig.Font_size = 8; 

SetFigure(cfg_fig, gcf,1)

%% convert to table data for anovas

 % group = [reshape(zeros(size(data.day_mean_ctrl)), 1, numel(data.day_mean_ctrl)), reshape(ones(size(data.day_mean_tau)), 1, numel(data.day_mean_tau))]; 
A_tbl = table(repmat([data.Jan25.sub' data.Jul24.sub'],1,5)',...
    repmat([data.Jan25.geno' data.Jul24.geno'],1,5)',...
    [zeros(1,length([data.Jan25.geno' data.Jul24.geno']))+1,zeros(1,length([data.Jan25.geno' data.Jul24.geno']))+2, zeros(1,length([data.Jan25.geno' data.Jul24.geno']))+3, zeros(1,length([data.Jan25.geno' data.Jul24.geno']))+4, zeros(1,length([data.Jan25.geno' data.Jul24.geno']))+5]',...
    [data.Jan25.day_mean(:,1)' data.Jul24.day_mean(:,1)' , data.Jan25.day_mean(:,2)' data.Jul24.day_mean(:,2)', data.Jan25.day_mean(:,3)' data.Jul24.day_mean(:,3)', data.Jan25.day_mean(:,4)' data.Jul24.day_mean(:,4)', data.Jan25.day_mean(:,5)' data.Jul24.day_mean(:,5)']',...
    [data.Jan25.day_first(:,1)' data.Jul24.day_first(:,1)' , data.Jan25.day_first(:,2)' data.Jul24.day_first(:,2)', data.Jan25.day_first(:,3)' data.Jul24.day_first(:,3)', data.Jan25.day_first(:,4)' data.Jul24.day_first(:,4)', data.Jan25.day_first(:,5)' data.Jul24.day_first(:,5)']',...
    'VariableNames',{'Subject', 'Geno','Day', 'Day_mean', 'Day_first'}); 

A_tbl.Day = categorical(A_tbl.Day); 
A_tbl.Geno = categorical(A_tbl.Geno); 

% stats

% r anova

% model
day_mdl = fitrm(A_tbl, 'Day_mean~Geno*Day', 'WithinDesign',table([1 2 3 4 5],'VariableNames', "Days"));
% repeated measures ANOVA
ranova_tbl = ranova(day_mdl, 'WithinModel', 'Days');
% results
disp(ranova_tbl);
%     multcompare(ranova_tbl)

% model
day_f_mdl = fitrm(A_tbl, 'Day_first~Geno*Day', 'WithinDesign',table([1 2 3 4 5],'VariableNames', "Days"));
% repeated measures ANOVA
ranova_df_tbl = ranova(day_f_mdl, 'WithinModel', 'Days');
% results
disp(ranova_df_tbl);

%% Collect the probe data

  figure(101)
    clf
%     title(coh{iC})
    subplot(1, 4, 1:3)
    cla
    hold on

    % for ii = 1:length(data.(coh{iC}).sub)
    %     plot(1:5, this_day_mean(ii,:))
    % 
    % end    
    plot(1:5, mean(data.day_mean_ctrl), 'color', c_ord(1,:), 'LineWidth',3)
    errorb(1:5, mean(data.day_mean_ctrl),MS_SEM_vec(data.day_mean_ctrl), 'color', c_ord(1,:), 'LineWidth',2)
    
    plot(1:5, mean(data.day_mean_tau), 'color', c_ord(2,:), 'LineWidth',3)
    errorb(1:5, mean(data.day_mean_tau),MS_SEM_vec(data.day_mean_tau), 'color', c_ord(2,:), 'LineWidth',2)


%     plot(1:5, mean(this_day_mean), 'k', 'LineWidth',2)
    % leg = data.(coh{iC}).sub; 
    % leg{end+1} = 'mean';
    legend({'Ctrl', 'Tau+'}, 'Box','off')
    ylabel('Mean escape latency')

    set(gca,'xtick', 1:5, 'XTickLabel', {'Day 1' 'Day2 ' 'Day 3' 'Day 4' 'Day 5'}, 'XTickLabelRotation', 45)

    subplot(1,4,4)
        
    b = MS_bar_w_err(data.cross_ctrl,data.cross_tau,c_ord(1:2,:), 1, 'ttest2' );
    b(1).FaceColor = 'none';
    b(1).EdgeColor = 'k';

    scatter(1+sort(MS_randn_range(length(data.cross_ctrl), 1, -.1, .1)), data.cross_ctrl,25,  c_ord(1,:), 'filled')
    scatter(2+sort(MS_randn_range(length(data.cross_tau), 1, -.1, .1)), data.cross_tau,25,  c_ord(2,:), 'filled')

    ylabel('Annulus crossings')
    set(gca, 'xtick', 1:2, 'xticklabel', {'Tau -', 'Tau +'}, 'XTickLabelRotation', 45)


SetFigure([], gcf,1)

%%

set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'Position',[50 50 1200 800]);
print(gcf, '-dpdf',[cd  filesep 'WM_summary.pdf'])
print(gcf, '-dpng',[cd  filesep 'WM_summary.png'])
