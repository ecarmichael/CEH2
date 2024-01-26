%% MASTER_Asmbly

if strcmp(computer, 'GLNXA64')
    
    codebase_dir = '/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared';
    ca_dir = '/home/williamslab/Documents/Github/CEH2';
    oasis_dir = '/home/williamslab/Documents/Github/OASIS_matlab';
    
    code_dir = '/home/williamslab/Documents/Github/Dos-Santos Assembly ICA/Dos-Santos Assembly ICA';
    
    RnR_dir = '/home/williamslab/Documents/Github/RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    data_dir = '/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter';
    
    
    
else
    
    codebase_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared';
    ca_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2';
    oasis_dir = 'C:\Users\ecarm\Documents\GitHub\OASIS_matlab';
    
    code_dir = 'C:\Users\ecarm\Downloads\Dos-Santos Assembly ICA\Dos-Santos Assembly ICA';
    
    RnR_dir = 'C:\Users\ecarm\Documents\GitHub\RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter';
    
end

restoredefaultpath
c_d = cd;
cd(oasis_dir)
addpath(genpath(oasis_dir));
oasis_setup

addpath(genpath(ca_dir));
addpath(genpath(codebase_dir))
addpath(genpath(RnR_dir));

addpath(code_dir)


cd(c_d)

%%

cd('C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter')
fig_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\checks';
% cd('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter')

f_list = dir('*data*');

move_thresh  = 8;
bin_size = [.25 .5];

out = [];
session = []; novel_idx = []; anx_idx = []; HS_idx = [];


for ii = 1:length(f_list)
    session{ii} = f_list(ii).name;
    
    % compute assemblies and related ReActs
    A_out{ii} = Pipeline_Asmbly(f_list(ii).name,bin_size, move_thresh);
    
    % Summary plots
    %     Pipline_Asmbly_plot(A_out{ii}, fig_dir);
    close all
    
    if ~isempty(strfind(f_list(ii).name, 'HATDS'))
        HS_idx(ii) = 1;
    else
        HS_idx(ii) = 0;
    end
    
    if ~isempty(strfind(f_list(ii).name, 'D1')) %|| ~isempty(strfind(f_list(ii).name, 'HATDS'))
        novel_idx(ii) = 1;
    end
    
    if ~isempty(strfind(f_list(ii).name, 'D5'))
        novel_idx(ii) = 0;
    end
    
    if ~isempty(strfind(f_list(ii).name, 'HAT'))
        anx_idx(ii) = 1;
    else
        anx_idx(ii) = 0;
    end
    
end

A_all = A_out;
A_out = A_out(1:11);


%%  Run again but for EVV data
data_dir = ('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV');
fig_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV\checks';
% cd('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter')
cd(data_dir);
j20_list = dir('*data*');



out = [];
J20_session = []; J20_novel_idx = []; D3_idx = [];
for ii = 1:length(j20_list)
    J20_session{ii} = j20_list(ii).name;
    cd(data_dir)
    
    J_out{ii} = Pipeline_Asmbly(j20_list(ii).name,bin_size, move_thresh);
    
    % Summary plots
    %     Pipline_Asmbly_plot(J_out{ii}, fig_dir);
    close all
    
    if ~isempty(strfind(j20_list(ii).name, 'D1')) %|| ~isempty(strfind(f_list(ii).name, 'HATDS'))
        J20_novel_idx(ii) = 1;
        D3_idx(ii) = 0;
    end
    
    if ~isempty(strfind(j20_list(ii).name, 'D3'))
        D3_idx(ii) = 1;
        J20_novel_idx(ii) = 0;
    end
    
    if ~isempty(strfind(j20_list(ii).name, 'D5'))
        J20_novel_idx(ii) = 0;
        D3_idx(ii) = 0;
    end
    
end

J_all = J_out;
A_out = J_out;
%% collect the data

Pre_n_Asmbly = []; Post_n_Asmbly = [];
Pre_r_Asmbly = []; Post_r_Asmbly = [];

S_Pre_n_Asmbly = []; S_Post_n_Asmbly = [];
S_Pre_r_Asmbly = []; S_Post_r_Asmbly = [];

for iB = length(bin_size):-1:1
    
    for iA = size(A_out,2):-1:1
        
        % pre
        Pre_n_Asmbly(iA,iB) = sum(A_out{iA}{iB}.REM_Pre_stats.p_val <0.05);
        
        Pre_r_Asmbly(iA,iB) = mean(A_out{iA}{iB}.REM_Pre_stats.rate(A_out{iA}{iB}.REM_Pre_stats.p_val <0.05));
        
        % post
        Post_n_Asmbly(iA,iB) = sum(A_out{iA}{iB}.REM_Post_stats.p_val <0.05);
        
        Post_r_Asmbly(iA,iB) = mean(A_out{iA}{iB}.REM_Post_stats.rate(A_out{iA}{iB}.REM_Post_stats.p_val <0.05));
        
        sub_list{iA} = A_out{iA}{iB}.info.subject;
        
        keep_idx = logical(A_out{iA}{iB}.REM_Pre_stats.p_val <0.05) & (A_out{iA}{iB}.REM_Post_stats.p_val <0.05);
        mean_ReAct_str(iA, iB) = mean(A_out{iA}{iB}.ReAct(keep_idx));
        
        
        % apply the same thing to assemblies with coherent spatial tuning
        % (only using zscore of the varience of the centroid position.
        
        
        
        % pre
        these_z_cent = []; these_z_peak = [];
        for ii = length(A_out{iA}{iB}.map):-1:1
            these_z_cent(ii) = A_out{iA}{iB}.map{ii}.cent_z;
            these_z_peak(ii) = A_out{iA}{iB}.map{ii}.peak_z;
        end
        these_sig = A_out{iA}{iB}.REM_Pre_stats.p_val;
        
        S_Pre_n_Asmbly(iA,iB) = sum((these_sig <0.05) & (these_z_cent < -1.96));
        S_Pre_r_Asmbly(iA,iB) = mean(A_out{iA}{iB}.REM_Pre_stats.rate((these_sig <0.05) & (these_z_cent < -1.96)));
        
        % post
        these_sig = A_out{iA}{iB}.REM_Post_stats.p_val;
        
        S_Post_n_Asmbly(iA,iB) = sum((these_sig <0.05) & (these_z_cent < -1.96));
        S_Post_r_Asmbly(iA,iB) = mean(A_out{iA}{iB}.REM_Post_stats.rate((these_sig <0.05) & (these_z_cent < -1.96)));
        
        
        keep_idx = logical(A_out{iA}{iB}.REM_Pre_stats.p_val <0.05) & (A_out{iA}{iB}.REM_Post_stats.p_val <0.05);
        mean_S_ReAct_str(iA, iB) = mean(A_out{iA}{iB}.ReAct(keep_idx & (these_z_cent < -1.96)));
        
    end
    
end

%% plot the number of Assemblies pre and post.

max_n_A = max([mean(Pre_n_Asmbly) ; mean(Post_n_Asmbly) ]);
max_n_A = max_n_A*1.5;
max_r_A = max([mean(Pre_r_Asmbly) ; mean(Post_r_Asmbly) ]);
max_r_A = max_r_A*1.5;

max_n_SA = max([mean(S_Pre_n_Asmbly) ; mean(S_Post_n_Asmbly) ]);
max_n_SA = max_n_SA*1.5;
max_r_SA = max([mean(S_Pre_r_Asmbly) ; mean(S_Post_r_Asmbly) ]);
max_r_SA = max_r_SA*1.5;


c_ord = MS_linspecer(5);
c_ord_s = c_ord*.7;

HS_idx = logical(HS_idx(1:length(A_out)));
n_idx = logical(novel_idx(1:length(A_out)));
a_idx = logical(anx_idx(1:length(A_out)));

n_idx = n_idx & ~HS_idx;
a_idx = a_idx & ~HS_idx;
f_idx = ~n_idx & ~HS_idx; 

lt1_idx = n_idx & ~a_idx & ~HS_idx;
lt5_idx = ~n_idx & ~a_idx & ~HS_idx;

H1_idx = n_idx & a_idx & ~HS_idx;
H5_idx = ~n_idx & a_idx & ~HS_idx;

% j20 =
n_idx = logical(novel_idx(1:length(A_out)));
a_idx = logical(novel_idx(1:length(A_out)));
n_idx = logical(J20_novel_idx);



for iB = length(bin_size):-1:1
    n = 4; m = 5;
    figure(100 + iB)
    
    subplot(n,m,1)
    cla
    MS_bar_w_err(Pre_n_Asmbly(n_idx,iB), Post_n_Asmbly(n_idx,iB),c_ord(1,:))
    
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title(['Novel (' num2str(bin_size(iB)) 's bins)'])
    ylabel('# assemblies')
    ylim([0 max_n_A(iB)])
    
    
    subplot(n,m,2)
    cla
    MS_bar_w_err(Pre_n_Asmbly(f_idx,iB), Post_n_Asmbly(f_idx,iB),c_ord(2,:))
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Fam')
    ylim([0 max_n_A(iB)])
    
    
    subplot(n,m,3)
    cla
    MS_bar_w_err(Pre_n_Asmbly(~a_idx,iB), Post_n_Asmbly(~a_idx,iB),c_ord(3,:))
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Linear')
    ylim([0 max_n_A(iB)])
    
    subplot(n,m,4)
    cla
    MS_bar_w_err(Pre_n_Asmbly(a_idx,iB), Post_n_Asmbly(a_idx,iB),c_ord(4,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Anxiety')
    ylim([0 max_n_A(iB)])
    
    subplot(n,m,5)
    cla
    MS_bar_w_err(Pre_n_Asmbly(HS_idx,iB), Post_n_Asmbly(HS_idx,iB),c_ord(5,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('A. Switch')
    ylim([0 max_n_A(iB)])
    
    
    % same for the rate
    
    subplot(n,m,m+1)
    cla
    MS_bar_w_err(Pre_r_Asmbly(n_idx,iB), Post_r_Asmbly(n_idx,iB),c_ord(1,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title(['Novel (' num2str(bin_size(iB)) 's bins)'])
    ylabel('ReAct/min')
    ylim([0 max_r_A(iB)])
    
    
    
    subplot(n,m,m+2)
    cla
    MS_bar_w_err(Pre_r_Asmbly(f_idx,iB), Post_r_Asmbly(f_idx,iB),c_ord(2,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Fam')
    ylim([0 max_r_A(iB)])
    
    
    subplot(n,m,m+3)
    cla
    MS_bar_w_err(Pre_r_Asmbly(~a_idx,iB), Post_r_Asmbly(~a_idx,iB),c_ord(3,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Linear')
    ylim([0 max_r_A(iB)])
    
    
    subplot(n,m,m+4)
    cla
    MS_bar_w_err(Pre_r_Asmbly(a_idx,iB), Post_r_Asmbly(a_idx,iB),c_ord(4,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Anxiety')
    ylim([0 max_r_A(iB)])
    
    
    subplot(n,m,m+5)
    cla
    MS_bar_w_err(Pre_r_Asmbly(HS_idx,iB), Post_r_Asmbly(HS_idx,iB),c_ord(5,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('A. Switch')
    ylim([0 max_r_A(iB)])
    
    
    
    % repeat for assemblies with coherent spatial tuning.
    subplot(n,m,(m*2)+1)
    cla
    MS_bar_w_err(S_Pre_n_Asmbly(n_idx,iB), S_Post_n_Asmbly(n_idx,iB),c_ord_s(1,:))
    
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title(['Novel (' num2str(bin_size(iB)) 's bins)'])
    ylabel('# assemblies')
    ylim([0 max_n_SA(iB)])
    
    
    subplot(n,m,(m*2)+2)
    cla
    MS_bar_w_err(S_Pre_n_Asmbly(f_idx,iB), S_Post_n_Asmbly(f_idx,iB),c_ord_s(2,:))
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Fam')
    ylim([0 max_n_SA(iB)])
    
    
    subplot(n,m,(m*2)+3)
    cla
    MS_bar_w_err(S_Pre_n_Asmbly(~a_idx,iB), S_Post_n_Asmbly(~a_idx,iB),c_ord_s(3,:))
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Linear')
    ylim([0 max_n_SA(iB)])
    
    subplot(n,m,(m*2)+4)
    cla
    MS_bar_w_err(S_Pre_n_Asmbly(a_idx,iB), S_Post_n_Asmbly(a_idx,iB),c_ord_s(4,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Anxiety')
    ylim([0 max_n_SA(iB)])
    
    subplot(n,m,(m*2)+5)
    cla
    MS_bar_w_err(S_Pre_n_Asmbly(HS_idx,iB), S_Post_n_Asmbly(HS_idx,iB),c_ord_s(5,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('A. Switch')
    ylim([0 max_n_SA(iB)])
    
    
    % same for the rate
    
    subplot(n,m,(m*3)+1)
    cla
    MS_bar_w_err(S_Pre_r_Asmbly(n_idx,iB), S_Post_r_Asmbly(n_idx,iB),c_ord_s(1,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title(['Novel (' num2str(bin_size(iB)) 's bins)'])
    ylabel('ReAct/min')
    ylim([0 max_r_SA(iB)])
    
    
    
    subplot(n,m,(m*3)+2)
    cla
    MS_bar_w_err(S_Pre_r_Asmbly(f_idx,iB), S_Post_r_Asmbly(f_idx,iB),c_ord_s(2,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Fam')
    ylim([0 max_r_SA(iB)])
    
    
    subplot(n,m,(m*3)+3)
    cla
    MS_bar_w_err(S_Pre_r_Asmbly(~a_idx,iB), S_Post_r_Asmbly(~a_idx,iB),c_ord_s(3,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Linear')
    ylim([0 max_r_SA(iB)])
    
    
    subplot(n,m,(m*3)+4)
    cla
    MS_bar_w_err(S_Pre_r_Asmbly(a_idx,iB), S_Post_r_Asmbly(a_idx,iB),c_ord_s(4,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Anxiety')
    ylim([0 max_r_SA(iB)])
    
    
    subplot(n,m,(m*3)+5)
    cla
    MS_bar_w_err(S_Pre_r_Asmbly(HS_idx,iB), S_Post_r_Asmbly(HS_idx,iB),c_ord_s(5,:));
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('A. Switch')
    ylim([0 max_r_SA(iB)])
    
    
end
%% collect j20 and control reactivation strength
A_out = A_all(1:11); 

A_Pre_n = []; A_Post_n = [];
A_Pre_r = []; A_Post_r = [];
A_cent = []; A_peak = []; 
A_ReAct_all = []; 

SA_Pre_n = []; SA_Post_n = [];
SA_Pre_r = []; SA_Post_r = [];
SA_cent = []; SA_peak = []; 

J_Pre_n = []; J_Post_n = [];
J_Pre_r = []; J_Post_r = [];
J_cent = []; J_peak = []; 
J_ReAct_all = []; 

SJ_Pre_n = []; SJ_Post_n = [];
SJ_Pre_r = []; SJ_Post_r = [];
SJ_cent = []; SJ_peak = []; 

A_ReAct = []; J_ReAct = []; 
SA_ReAct = []; SJ_ReAct = []; 

A_sub_list = []; J_sub_list = [];


for iB = length(bin_size):-1:1
%     all_cent = 
    for iA = size(A_out,2):-1:1
        
        % pre
        A_Pre_n(iA,iB) = sum(A_out{iA}{iB}.REM_Pre_stats.p_val <0.05);
        
        A_Pre_r(iA,iB) = mean(A_out{iA}{iB}.REM_Pre_stats.rate(A_out{iA}{iB}.REM_Pre_stats.p_val <0.05));
        
        % post
        A_Post_a(iA,iB) = sum(A_out{iA}{iB}.REM_Post_stats.p_val <0.05);
        
        A_Post_r(iA,iB) = mean(A_out{iA}{iB}.REM_Post_stats.rate(A_out{iA}{iB}.REM_Post_stats.p_val <0.05));
        
        A_sub_list{iA} = A_out{iA}{iB}.info.subject;
        
        keep_idx = logical(A_out{iA}{iB}.REM_Pre_stats.p_val <0.05) & (A_out{iA}{iB}.REM_Post_stats.p_val <0.05);
        A_ReAct(iA, iB) = mean(A_out{iA}{iB}.ReAct(keep_idx));
        
                A_ReAct_all{iB}{iA} = mean(A_out{iA}{iB}.REM_Post_proj,2) - mean(A_out{iA}{iB}.REM_Pre_proj,2); 

% Coherent spatial maps
these_z_cent = []; these_z_peak = [];
        for ii = length(A_out{iA}{iB}.map):-1:1
            these_z_cent(ii) = A_out{iA}{iB}.map{ii}.cent_z;
            these_z_peak(ii) = A_out{iA}{iB}.map{ii}.peak_z;
            A_cent{iB}{iA}(ii) = A_out{iA}{iB}.map{ii}.cent_z;
            A_peak{iB}{iA}(ii) = A_out{iA}{iB}.map{ii}.peak_z;
        end
        

        these_sig = A_out{iA}{iB}.REM_Pre_stats.p_val;
        
        SA_Pre_n(iA,iB) = sum((these_sig <0.05) & (these_z_cent < -1.96));
        SA_Pre_r(iA,iB) = mean(A_out{iA}{iB}.REM_Pre_stats.rate((these_sig <0.05) & (these_z_cent < -1.96)));
        
        % post
        these_sig = A_out{iA}{iB}.REM_Post_stats.p_val;
        
        SA_Post_n(iA,iB) = sum((these_sig <0.05) & (these_z_cent < -1.96));
        SA_Post_r(iA,iB) = mean(A_out{iA}{iB}.REM_Post_stats.rate((these_sig <0.05) & (these_z_cent < -1.96)));
        
        
        keep_idx = logical(A_out{iA}{iB}.REM_Pre_stats.p_val <0.05) & (A_out{iA}{iB}.REM_Post_stats.p_val <0.05);
        SA_ReAct(iA, iB) = mean(A_out{iA}{iB}.ReAct(keep_idx & (these_z_cent < -1.96)));
    end
end


for iB = length(bin_size):-1:1
    for iA = size(J_out,2):-1:1
        
        % pre
        J_Pre_n(iA,iB) = sum(J_out{iA}{iB}.REM_Pre_stats.p_val <0.05);
        
        J_Pre_r(iA,iB) = mean(J_out{iA}{iB}.REM_Pre_stats.rate(J_out{iA}{iB}.REM_Pre_stats.p_val <0.05));
        
        % post
        J_Post_a(iA,iB) = sum(J_out{iA}{iB}.REM_Post_stats.p_val <0.05);
        
        J_Post_r(iA,iB) = mean(J_out{iA}{iB}.REM_Post_stats.rate(J_out{iA}{iB}.REM_Post_stats.p_val <0.05));
        
        J_sub_list{iA} = J_out{iA}{iB}.info.subject;
        
        keep_idx = logical(J_out{iA}{iB}.REM_Pre_stats.p_val <0.05) & (J_out{iA}{iB}.REM_Post_stats.p_val <0.05);
        J_ReAct(iA, iB) = mean(J_out{iA}{iB}.ReAct(keep_idx));
        
        J_ReAct_all{iB}{iA} = mean(J_out{iA}{iB}.REM_Post_proj,2) - mean(J_out{iA}{iB}.REM_Pre_proj,2); 
        
        % Coherent spatial maps
        these_z_cent = []; these_z_peak = [];
        for ii = length(J_out{iA}{iB}.map):-1:1
            these_z_cent(ii) = J_out{iA}{iB}.map{ii}.cent_z;
            these_z_peak(ii) = J_out{iA}{iB}.map{ii}.peak_z;
            J_cent{iB}{iA}(ii) = J_out{iA}{iB}.map{ii}.cent_z;
            J_peak{iB}{iA}(ii) = J_out{iA}{iB}.map{ii}.peak_z;
        end
        these_sig = J_out{iA}{iB}.REM_Pre_stats.p_val;
        
        SJ_Pre_n(iA,iB) = sum((these_sig <0.05) & (these_z_cent < -1.96));
        SJ_Pre_r(iA,iB) = mean(J_out{iA}{iB}.REM_Pre_stats.rate((these_sig <0.05) & (these_z_cent < -1.96)));
        
        % post
        these_sig = J_out{iA}{iB}.REM_Post_stats.p_val;
        
        SJ_Post_n(iA,iB) = sum((these_sig <0.05) & (these_z_cent < -1.96));
        SJ_Post_r(iA,iB) = mean(J_out{iA}{iB}.REM_Post_stats.rate((these_sig <0.05) & (these_z_cent < -1.96)));
        
        
        keep_idx = logical(J_out{iA}{iB}.REM_Pre_stats.p_val <0.05) & (J_out{iA}{iB}.REM_Post_stats.p_val <0.05);
        SJ_ReAct(iA, iB) = mean(J_out{iA}{iB}.ReAct(keep_idx & (these_z_cent < -1.96)));
        
%         SJ_cent{iB}{iA} = J_cent{iB}(keep_idx); 
        
    end
end
%% plot the J20 vs control novel VS familiar 
A_out = A_all(1:11); 

HS_idx = logical(HS_idx(1:length(A_out)));
n_idx = logical(novel_idx(1:length(A_out)));
a_idx = logical(anx_idx(1:length(A_out)));

J_n_idx = logical(J20_novel_idx); 

figure(4000)


hb = bar([nanmean(A_ReAct(n_idx, iB)), nanmean(J_ReAct(J_n_idx, iB))]', 'FaceColor', c_ord(1,:), 'EdgeColor', c_ord(1,:));
hold on
eb = errorbar([nanmean(data_a), nanmean(data_b)], [MS_SEM(data_a) ,MS_SEM(data_b)]);
eb.LineStyle = 'none';
eb.Color = 'k';


%% check the assmbly reactivation strength against centorid Z
% control
HS_idx = logical(HS_idx(1:length(A_out)));
n_idx = logical(novel_idx(1:length(A_out)));
a_idx = logical(anx_idx(1:length(A_out)));

n_idx = n_idx & ~HS_idx;
a_idx = a_idx & ~HS_idx;
f_idx = ~n_idx & ~HS_idx; 

lt1_idx = n_idx & ~a_idx & ~HS_idx;
lt5_idx = ~n_idx & ~a_idx & ~HS_idx;

H1_idx = n_idx & a_idx & ~HS_idx;
H5_idx = ~n_idx & a_idx & ~HS_idx;

A_n_cent = []; A_n_peak = []; A_n_ReAct = []; A_n_sub = [];
A_f_cent = []; A_f_peak = []; A_f_ReAct = []; A_f_sub = [];

A_h_cent = []; A_h_peak = []; A_h_ReAct = []; A_h_sub = [];
A_l_cent = []; A_l_peak = []; A_l_ReAct = []; A_l_sub = [];

for ii = 1:length(A_out)
    
    if ismember(ii, find(n_idx))
        A_n_cent = [A_n_cent A_cent{iB}{ii}];
        A_n_peak = [A_n_peak A_peak{iB}{ii}];
        A_n_ReAct = [A_n_ReAct A_ReAct_all{iB}{ii}'];
        a_n_sub = [A_n_sub repmat(ii, 1, length(A_ReAct_all{iB}{ii}))];
    end
    
    if ismember(ii, find(f_idx))
        A_f_cent = [A_f_cent A_cent{iB}{ii}];
        A_f_peak = [A_f_peak A_peak{iB}{ii}];
        A_f_ReAct = [A_f_ReAct A_ReAct_all{iB}{ii}'];
        A_f_sub = [A_f_sub repmat(ii, 1, length(A_ReAct_all{iB}{ii}))];
    end
    
end

%% J20
D5_idx = ~J20_novel_idx & ~D3_idx;
J_n_cent = []; J_n_peak = []; J_n_ReAct = []; J_n_sub = [];
J_f_cent = []; J_f_peak = []; J_f_ReAct = []; J_f_sub = [];

for ii = 1:length(J_out)
    
    if ismember(ii, find(J20_novel_idx))
        J_n_cent = [J_n_cent J_cent{iB}{ii}];
        J_n_peak = [J_n_peak J_peak{iB}{ii}];
        J_n_ReAct = [J_n_ReAct J_ReAct_all{iB}{ii}'];
        J_n_sub = [J_n_sub repmat(ii, 1, length(J_ReAct_all{iB}{ii}))];
    end
    
    if ismember(ii, find(D5_idx))
        J_f_cent = [J_f_cent J_cent{iB}{ii}];
        J_f_peak = [J_f_peak J_peak{iB}{ii}];
        J_f_ReAct = [J_f_ReAct J_ReAct_all{iB}{ii}'];
        J_f_sub = [J_f_sub repmat(ii, 1, length(J_ReAct_all{iB}{ii}))];
    end
    
end
%%
 % control
figure(5000)
subplot(2,2,1)
scatter(A_n_cent, A_n_ReAct, 35, c_ord(1,:), 'filled');
[c, p] = corrcoef(A_n_cent, A_n_ReAct); 
text(-1, .7, ['r = ' num2str(c(1,2),2) ' p= ' num2str(p(1,2),2)], 'FontWeight', 'bold')
xlabel('(<- coherent) Z Map centroid')
ylabel('ReAct Str.')
title('Control novel')
ylim([-1 1]);
xlim([-5 5]);
lsline

subplot(2,2,2)
scatter(A_f_cent, A_f_ReAct, 35, c_ord(2,:), 'filled');
[c, p] = corrcoef(A_f_cent, A_f_ReAct); 
text(-1, .7, ['r = ' num2str(c(1,2),2) ' p= ' num2str(p(1,2),2)], 'FontWeight', 'bold')
xlabel('(<- coherent) Z Map centroid')
ylabel('ReAct Str.')
title('familiar')
ylim([-1 1]);
xlim([-5 5]);
lsline

subplot(2,2,3)
scatter(A_n_peak, A_n_ReAct, 35, c_ord(1,:), 'filled');
[c, p] = corrcoef(A_n_peak, A_n_ReAct); 
text(-1, .7, ['r = ' num2str(c(1,2),2) ' p= ' num2str(p(1,2),2)], 'FontWeight', 'bold')
ylabel('ReAct Str.')
xlabel('(<- coherent) Z Map peak')
ylim([-1 1]);
xlim([-5 5]);
lsline

subplot(2,2,4)
scatter(A_f_peak, A_f_ReAct, 35, c_ord(2,:), 'filled');
[c, p] = corrcoef(A_f_peak, A_f_ReAct); 
text(-1, .7, ['r = ' num2str(c(1,2),2) ' p= ' num2str(p(1,2),2)], 'FontWeight', 'bold')
xlabel('(<- coherent) Z Map peak')
ylabel('ReAct Str.')
ylim([-1 1]);
xlim([-5 5]);
lsline

% J20
figure(5001)
subplot(2,2,1)
scatter(J_n_cent, J_n_ReAct, 35, c_ord_s(1,:), 'filled');
[c, p] = corrcoef(J_n_cent, J_n_ReAct); 
text(-1, .7, ['r = ' num2str(c(1,2),2) ' p= ' num2str(p(1,2),2)], 'FontWeight', 'bold')
xlabel('(<- coherent) Z Map centroid')
ylabel('ReAct Str.')
title('J20 novel')
ylim([-1 1]);
xlim([-5 5]);
lsline

subplot(2,2,2)
scatter(J_f_cent, J_f_ReAct, 35, c_ord_s(2,:), 'filled');
[c, p] = corrcoef(J_f_cent, J_f_ReAct); 
text(-1, .7, ['r = ' num2str(c(1,2),2) ' p= ' num2str(p(1,2),2)], 'FontWeight', 'bold')
xlabel('(<- coherent) Z Map centroid')
ylabel('ReAct Str.')
title('familiar')
ylim([-1 1]);
xlim([-5 5]);
lsline

subplot(2,2,3)
scatter(J_n_peak, J_n_ReAct, 35, c_ord_s(1,:), 'filled');
[c, p] = corrcoef(J_n_peak, J_n_ReAct); 
text(-1, .7, ['r = ' num2str(c(1,2),2) ' p= ' num2str(p(1,2),2)], 'FontWeight', 'bold')
ylabel('ReAct Str.')
xlabel('(<- coherent) Z Map peak')
ylim([-1 1]);
xlim([-5 5]);
lsline

subplot(2,2,4)
ax = scatter(J_f_peak, J_f_ReAct, 35, c_ord_s(2,:), 'filled');
lsline
[c, p] = corrcoef(J_f_peak, J_f_ReAct); 
text(-1, .7, ['r = ' num2str(c(1,2),2) ' p= ' num2str(p(1,2),2)], 'FontWeight', 'bold')
xlabel('(<- coherent) Z Map peak')
ylabel('ReAct Str.')
ylim([-1 1]);
xlim([-5 5]);

%% to do

% work out the location of the reactivated assemblies. Is there a bias? 
