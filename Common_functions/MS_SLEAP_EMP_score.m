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
    f_out.(f_list(iF).name(1:end-4)).metrics.CT = sum(out.emp_idx(1:end-1) == 1 & out.emp_idx(2:end) == 3); % count closed to transition
    f_out.(f_list(iF).name(1:end-4)).metrics.T_O =sum(out.emp_idx(1:end-1) == 3 & out.emp_idx(2:end) == 2); % count closed to transition
    
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
%% simple plots


figure(909)
subplot(2,1,1)
plot(1:length(table_out.Open_prct), table_out.Open_prct, '.b', 'markersize', 20)
set(gca, 'XTick', 1:length(table_out.Open_prct), 'XTickLabel', table_out.Subject, 'XTickLabelRotation', 65)
ylabel('prct time in Open arm')

subplot(2,1,2)
plot(1:length(table_out.Open_prct), table_out.T2O, '.r', 'markersize', 20)
set(gca, 'XTick', 1:length(table_out.Open_prct), 'XTickLabel', table_out.Subject, 'XTickLabelRotation', 65)
ylabel('Number of Open arm entries')


%% split based on aTau
table.geno = logical([0     0     0     0     0     1     0     1     1     0     0     0     1     1     0     1     1     0     0     0     0     0     0     0     0     0     0]); 

c_ord = MS_linspecer(4);
figure(909)
clf
hold on
scatter(1+sort(MS_randn_range(length(table_out.Closed_prct(table_out.Geno == 0)), 1, -.1, .1)), table_out.Closed_prct(table_out.Geno == 0),25,  c_ord(1,:), 'filled')
scatter(2+sort(MS_randn_range(length(table_out.Closed_prct(table_out.Geno == 1)), 1, -.1, .1)), table_out.Closed_prct(table_out.Geno == 1),25,  c_ord(2,:), 'filled')

scatter(5+sort(MS_randn_range(length(table_out.T2O(table_out.Geno == 0)), 1, -.1, .1)), table_out.T2O(table_out.Geno == 0),25,  c_ord(1,:), 'filled')
scatter(6+sort(MS_randn_range(length(table_out.T2O(table_out.Geno == 1)), 1, -.1, .1)), table_out.T2O(table_out.Geno == 1),25,  c_ord(2,:), 'filled')

[hb h, p]= MS_bar_w_err(table_out.Closed_prct(table_out.Geno == 0), table_out.Closed_prct(table_out.Geno == 1),c_ord(1,:),  1, 'ttest2', 1:2)

hb(1).FaceColor = 'none';
hb(1).EdgeColor = 'k';

[hb h, p]= MS_bar_w_err(table_out.T2O(table_out.Geno == 0), table_out.T2O(table_out.Geno == 1),c_ord(2,:),  1, 'ttest2', 5:6)
hb(1).FaceColor = 'none';
hb(1).EdgeColor = 'k';

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