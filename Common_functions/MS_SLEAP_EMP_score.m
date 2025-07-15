function [f_out] = MS_SLEAP_EMP_score(data_dir)
%% MS_SLEAP_EMP_score: collects the state indicies from MS_SLEAP_EMP and converts it to a csv with all the relevant metrics. 
%
%
%
%    Inputs: 
%    - data_dir: [path] path to intermediate files. 
%
%
%
%    Outputs: 
%    - f_out: [file containing a table of the EMP metrics]
%
%
%
%
% EC 2024-12-19   initial version 
%
%
%
%% initialize

if nargin < 1
    data_dir = cd;
end


%% collect all the data;
cd(data_dir)

f_list = dir('*.mat'); 

for iF = 1:length(f_list)
    
    load(f_list(iF).name); 
    
%     
%     c_idx = find(contains(out.labels, 'closed')); 
%     c_idx = find(contains(out.labels, 'closed')); 
%     c_idx = find(contains(out.labels, 'closed')); 
%     c_idx = find(contains(out.labels, 'closed')); 

    
        
    f_out.(f_list(iF).name(1:end-4)).metrics.closed = (sum(out.emp_idx == 1)./length(out.emp_idx))*100; 
    f_out.(f_list(iF).name(1:end-4)).metrics.open = (sum(out.emp_idx == 2 | out.emp_idx == 4)./length(out.emp_idx))*100; 
    f_out.(f_list(iF).name(1:end-4)).metrics.trans = (sum(out.emp_idx == 3)./length(out.emp_idx))*100; 

    
    f_out.(f_list(iF).name(1:end-4)).metrics.N_head_dip = length(out.dip_IV.tstart); 
    
    % transitions. 
    f_out.(f_list(iF).name(1:end-4)).metrics.CT = length(out.CtT_IV.tstart); % sum(out.emp_idx(1:end-1) == 1 & out.emp_idx(2:end) == 3); % count closed to transition
    f_out.(f_list(iF).name(1:end-4)).metrics.T_O =length(out.OtT_IV.tstart); % sum(out.emp_idx(1:end-1) == 3 & out.emp_idx(2:end) == 2); % count closed to transition
    
    
    
    % compile for the table; 
    subject{iF} =  strrep(f_list(iF).name(1:end-15), '_', ' ');
    date_id{iF} = f_list(iF).name(end-13:end-4); 
    Closed_prct(iF) = f_out.(f_list(iF).name(1:end-4)).metrics.closed;
    Open_prct(iF) = f_out.(f_list(iF).name(1:end-4)).metrics.open;
    Trans_prct(iF) = f_out.(f_list(iF).name(1:end-4)).metrics.trans;
    N_head_dip(iF) = f_out.(f_list(iF).name(1:end-4)).metrics.N_head_dip;
    C2T(iF) = f_out.(f_list(iF).name(1:end-4)).metrics.CT;
    T2O(iF) = f_out.(f_list(iF).name(1:end-4)).metrics.T_O;


end

% convert to table
table_out = table(subject', date_id', Closed_prct', Open_prct', Trans_prct', N_head_dip', C2T', T2O','VariableNames',  ["Subject", "Date",  "Closed_prct", "Open_prct", "Trans_prct", "N_head_dip", "C2T", "T2O"]);



%% split based on aTau
table_out.geno = logical([0     0     0     0     0     1     0     1     1     0     0     0     1     1     0, ...
    1     1     0     0     0     0     0     0     0     0     0     ])'; 

c_ord = MS_linspecer(4);
figure(909)
clf
title('EPM')
subplot(2,2,1)
hold on

[h_c, p, stats] = MS_rain_plot(table_out.Closed_prct, table_out.geno, c_ord(1:2,:), 'ttest2');
fprintf('Closed ttest: df: %2.0f  tstat: %2.2f p: %2.3f\n', stats.df, stats.tstat, p)
% scatter(1+sort(MS_randn_range(length(table_out.Closed_prct(table_out.geno == 0)), 1, -.1, .1)), table_out.Closed_prct(table_out.geno == 0),25,  c_ord(1,:), 'filled')
% scatter(2+sort(MS_randn_range(length(table_out.Closed_prct(table_out.geno == 1)), 1, -.1, .1)), table_out.Closed_prct(table_out.geno == 1),25,  c_ord(2,:), 'filled')

% scatter(4+sort(MS_randn_range(length(table_out.T2O(table_out.geno == 0)), 1, -.1, .1)), table_out.T2O(table_out.geno == 0),25,  c_ord(1,:), 'filled')
% scatter(5+sort(MS_randn_range(length(table_out.T2O(table_out.geno == 1)), 1, -.1, .1)), table_out.T2O(table_out.geno == 1),25,  c_ord(2,:), 'filled')

% [hb h, p]= MS_bar_w_err(table_out.Closed_prct(table_out.geno == 0), table_out.Closed_prct(table_out.geno == 1),c_ord(1,:),  1, 'ttest2', 1:2)

% [hb, h, p]
% 
% hb(1).FaceColor = 'none';
% hb(1).EdgeColor = 'k';

% [hb h, p]= MS_bar_w_err(table_out.T2O(table_out.geno == 0), table_out.T2O(table_out.geno == 1),c_ord(2,:),  1, 'ttest2', 4:5)

xlabel('% of time in closed arm')

set(gca, 'ytick', 1:2, 'yticklabel', {'Tau -', 'Tau +'}, 'YTickLabelRotation' , 45)
ylim([0 3])
xlim([0 100])


subplot(2,2,3)
hold on
[h_c, p, stats] = MS_rain_plot(table_out.C2T, table_out.geno, c_ord(1:2,:), 'ttest2');
fprintf('Closed-to-open transitions ttest: df: %2.0f  tstat: %2.2f p: %2.3f\n', stats.df, stats.tstat, p)

xlabel('Number of transition entries')
set(gca, 'ytick', 1:2, 'yticklabel', {'Tau -', 'Tau +'}, 'YTickLabelRotation' , 45)
% boxplot(table_out.Closed_prct,table_out.geno)
ylim([0 3])
x_lim = xlim; 

xlim([0 x_lim(2)])


subplot(2,2,2)
[h_c, p, stats] = MS_rain_plot(table_out.Open_prct, table_out.geno, c_ord(1:2,:), 'ttest2');
fprintf('Open %% ttest: df: %2.0f  tstat: %2.2f p: %2.3f\n', stats.df, stats.tstat, p)

xlabel('% of time in open arm')
set(gca, 'ytick', 1:2, 'yticklabel', {'Tau -', 'Tau +'}, 'YTickLabelRotation' , 45)
% x = 1:length(table_out.Open_prct);
% plot(x(table_out.geno == 0), table_out.Open_prct(table_out.geno == 0), '.b', 'markersize', 40, 'filled')
% set(gca, 'XTick', 1:length(table_out.Open_prct), 'XTickLabel', table_out.Subject, 'XTickLabelRotation', 65)
% ylabel('prct time in Open arm')
ylim([0 3])
xlim([0 100])

subplot(2,2,4)
[h_c, p, stats] = MS_rain_plot(table_out.T2O, table_out.geno, c_ord(1:2,:), 'ttest2');
fprintf('transition-to-open transitions ttest: df: %2.0f  tstat: %2.2f p: %2.3f\n', stats.df, stats.tstat, p)

xlabel('Number of open arm entries')
set(gca, 'ytick', 1:2, 'yticklabel', {'Tau -', 'Tau +'}, 'YTickLabelRotation' , 45)
% plot(1:length(table_out.Open_prct), table_out.T2O, '.r', 'markersize', 20)
% set(gca, 'XTick', 1:length(table_out.Open_prct), 'XTickLabel', table_out.Subject, 'XTickLabelRotation', 65)
% ylabel('Number of Open arm entries')
ylim([0 3])
x_lim = xlim; 
xlim([0 x_lim(2)])

SetFigure([], gcf, 1)
saveas(gcf, 'Summary_rain.png')
print(gcf, '-dpdf','Summary_rain.pdf')


%% sample plot but bars
table_out.geno = logical([0     0     0     0     0     1     0     1     1     0     0     0     1     1     0, ...
    1     1     0     0     0     0     0     0     0     0     0     ])'; 

c_ord = MS_linspecer(4);
figure(909)
clf
title('EPM')
subplot(2,2,1)
hold on

% [h_c, p, stats] = MS_rain_plot(table_out.Closed_prct, table_out.geno, c_ord(1:2,:), 'ttest2');
[h_c, eb, sc, p, stats] = MS_bar_w_err(table_out.Closed_prct(table_out.geno ==0),table_out.Closed_prct(table_out.geno ==1) , c_ord(1:2,:),1, 'ttest2', [1 2]);
fprintf('Closed ttest: df: %2.0f  tstat: %2.2f p: %2.3f\n', stats.df, stats.tstat, p)

h_c.FaceColor = 'none'; 
h_c.EdgeColor = 'k'; h_c.FaceColor = 'none'; 
h_c.EdgeColor = 'k'; 
% scatter(1+sort(MS_randn_range(length(table_out.Closed_prct(table_out.geno == 0)), 1, -.1, .1)), table_out.Closed_prct(table_out.geno == 0),25,  c_ord(1,:), 'filled')
% scatter(2+sort(MS_randn_range(length(table_out.Closed_prct(table_out.geno == 1)), 1, -.1, .1)), table_out.Closed_prct(table_out.geno == 1),25,  c_ord(2,:), 'filled')

% scatter(4+sort(MS_randn_range(length(table_out.T2O(table_out.geno == 0)), 1, -.1, .1)), table_out.T2O(table_out.geno == 0),25,  c_ord(1,:), 'filled')
% scatter(5+sort(MS_randn_range(length(table_out.T2O(table_out.geno == 1)), 1, -.1, .1)), table_out.T2O(table_out.geno == 1),25,  c_ord(2,:), 'filled')

% [hb h, p]= MS_bar_w_err(table_out.Closed_prct(table_out.geno == 0), table_out.Closed_prct(table_out.geno == 1),c_ord(1,:),  1, 'ttest2', 1:2)

% [hb, h, p]
% 
% hb(1).FaceColor = 'none';
% hb(1).EdgeColor = 'k';

% [hb h, p]= MS_bar_w_err(table_out.T2O(table_out.geno == 0), table_out.T2O(table_out.geno == 1),c_ord(2,:),  1, 'ttest2', 4:5)

ylabel('% of time in closed arm')

set(gca, 'xtick', 1:2, 'xticklabel', {'Tau -', 'Tau +'}, 'xTickLabelRotation' , 45)
xlim([-1 4])
y_lim = ylim; 

ylim([0 y_lim(2)])

subplot(2,2,3)
cla
hold on
[h_c, eb, sc, p, stats] = MS_bar_w_err(table_out.C2T(table_out.geno ==0),table_out.C2T(table_out.geno ==1) , c_ord(1:2,:),1, 'ttest2', [1 2]);
fprintf('Closed-to-open transitions ttest: df: %2.0f  tstat: %2.2f p: %2.3f\n', stats.df, stats.tstat, p)

h_c.FaceColor = 'none'; 
h_c.EdgeColor = 'k'; 

ylabel('Number of transition entries')
set(gca, 'xtick', 1:2, 'xticklabel', {'Tau -', 'Tau +'}, 'xTickLabelRotation' , 45)
% boxplot(table_out.Closed_prct,table_out.geno)
xlim([-1 4])
y_lim = ylim; 

ylim([0 y_lim(2)])


subplot(2,2,2)
cla
% [h_c, p, stats] = MS_rain_plot(table_out.Open_prct, table_out.geno, c_ord(1:2,:), 'ttest2');
[h_c, eb, sc, p, stats] = MS_bar_w_err(table_out.Open_prct(table_out.geno ==0),table_out.Open_prct(table_out.geno ==1) , c_ord(1:2,:),1, 'ttest2', [1 2]);
fprintf('Open %% ttest: df: %2.0f  tstat: %2.2f p: %2.3f\n', stats.df, stats.tstat, p)
h_c.FaceColor = 'none'; 
h_c.EdgeColor = 'k'; 

ylabel('% of time in open arm')
set(gca, 'xtick', 1:2, 'xticklabel', {'Tau -', 'Tau +'}, 'XTickLabelRotation' , 45)
% x = 1:length(table_out.Open_prct);
% plot(x(table_out.geno == 0), table_out.Open_prct(table_out.geno == 0), '.b', 'markersize', 40, 'filled')
% set(gca, 'XTick', 1:length(table_out.Open_prct), 'XTickLabel', table_out.Subject, 'XTickLabelRotation', 65)
% ylabel('prct time in Open arm')
xlim([-1 4])
y_lim = ylim; 

ylim([0 y_lim(2)])
subplot(2,2,4)
% [h_c, p, stats] = MS_rain_plot(table_out.T2O, table_out.geno, c_ord(1:2,:), 'ttest2');
[h_c, eb, sc, p, stats] = MS_bar_w_err(table_out.T2O(table_out.geno ==0),table_out.T2O(table_out.geno ==1) , c_ord(1:2,:),1, 'ttest2', [1 2]);
fprintf('transition-to-open transitions ttest: df: %2.0f  tstat: %2.2f p: %2.3f\n', stats.df, stats.tstat, p)
h_c.FaceColor = 'none'; 
h_c.EdgeColor = 'k'; 

ylabel('Number of open arm entries')
set(gca, 'xtick', 1:2, 'xticklabel', {'Tau -', 'Tau +'}, 'xTickLabelRotation' , 45)
% plot(1:length(table_out.Open_prct), table_out.T2O, '.r', 'markersize', 20)
% set(gca, 'XTick', 1:length(table_out.Open_prct), 'XTickLabel', table_out.Subject, 'XTickLabelRotation', 65)
% ylabel('Number of Open arm entries')
xlim([-1 4])
y_lim = ylim; 
ylim([0 y_lim(2)])

SetFigure([], gcf, 1)
saveas(gcf, 'Summary.png')
print(gcf, '-dpdf','Summary.pdf')

%% one shot to update the EPM data

cd(data_dir)
d = dir(data_dir);
[p_d] =fileparts(d(3).folder);


f_list = dir('*.mat'); 

for iF = 1:length(f_list)
    
    cd(data_dir);
    
    load(f_list(iF).name);

    
    s_name = f_list(iF).name(1:end-15);
%     d_name = f_list(iF).name(end-13:end-4); 
   h5_list = dir(fullfile([p_d filesep s_name],'**', '*.h5'));
   this_dir = h5_list(end).folder; 
   
   cd(this_dir); 
   
   MS_SLEAP_EMP(this_dir, save_dir, out.boxes)
   
    
    pause(2)
    close all
end

%% update metrics using precomputed data and boxes
cd(data_dir)
close all
f_list = dir('*.mat'); 

for iF = 1:length(f_list)
    
    cd(data_dir)

   load(f_list(iF).name)
   
   
      s_name = f_list(iF).name(1:end-15);
%     d_name = f_list(iF).name(end-13:end-4); 
   h5_list = dir(fullfile([p_d filesep s_name],'**', '*.h5'));
   this_dir = h5_list(end).folder; 
   
   cd(this_dir); 
   if ~isempty(dir('*DLC*'))
        MS_DLC_EMP(this_dir, 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\PoxR1\EPM\inter_temp', out.boxes);
   else
          MS_SLEAP_EMP(this_dir, 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\PoxR1\EPM\inter_temp', out.boxes);
   end
   pause(1)
   close all
   
    
end