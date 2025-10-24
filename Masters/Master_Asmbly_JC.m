%% MASTER_Asmbly

if strcmp(computer, 'GLNXA64')

    usr = char(java.lang.System.getProperty('user.name')); 
    
    codebase_dir = ['/home/' usr '/Documents/Github/vandermeerlab/code-matlab/shared'];
    ca_dir = ['/home/' usr '/Documents/Github/CEH2'];
    oasis_dir = ['/home/' usr '/Documents/Github/OASIS_matlab'];
    
    code_dir = ['/home/' usr '/Documents/Github/Dos-Santos Assembly ICA/Dos-Santos Assembly ICA'];
    
    RnR_dir = ['/home/' usr '/Documents/Github/RnR_methods'];
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    main_dir = ['/home/' usr '/'];
    
    
elseif strcmp(computer, 'MACA64')
    
    codebase_dir = '/Users/ecar/Documents/Github/vandermeerlab/code-matlab/shared';
    ca_dir = '/Users/ecar/Documents/Github/CEH2';
    oasis_dir = '/Users/ecar/Documents/Github/OASIS_matlab';
    
    code_dir = '/Users/ecar/Documents/Github/Dos-Santos Assembly ICA/Dos-Santos Assembly ICA';
    
    RnR_dir = '//Users/ecar/Documents/Github/RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    main_dir = '/Users/ecar/';
    
else
    
    codebase_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared';
    ca_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2';
    oasis_dir = 'C:\Users\ecarm\Documents\GitHub\OASIS_matlab';
    
    code_dir = 'C:\Users\ecarm\Downloads\Dos-Santos Assembly ICA\Dos-Santos Assembly ICA';
    
    RnR_dir = 'C:\Users\ecarm\Documents\GitHub\RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    main_dir = 'C:\Users\ecarm\';
    
end

restoredefaultpath
c_d = cd;
% cd(oasis_dir)
% addpath(genpath(oasis_dir));
% oasis_setup

addpath(genpath(ca_dir));
addpath(genpath(codebase_dir))
addpath(genpath(RnR_dir));

addpath(genpath(code_dir))


cd(c_d)


move_thresh  = 9;
bin_size = .5;

inter_dir = strrep([main_dir 'Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\PC9'], '\', filesep);
fig_dir = [inter_dir filesep 'checks'];

%%  Extract assembly data

% cd('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter')

cd(inter_dir);

f_list = dir('*data*');

A_out = []; P_out = [];
session = []; novel_idx = []; anx_idx = []; HS_idx = [];
method = 'binary';

for ii = 1:length(f_list)
    session{ii} = f_list(ii).name;
    % compute assemblies and related ReActs
    % A_out{ii} = Pipeline_Asmbly(f_list(ii).name,bin_size, move_thresh, method);
    % P_out{ii} = Pipeline_Asmbly_place(f_list(ii).name,bin_size, move_thresh, method);

    B_out{ii} = Pipeline_Asmbly_top_cells(f_list(ii).name,bin_size, move_thresh, method, [], 256);

    B_out{ii} = Pipeline_Asmbly_append_SWS(f_list(ii).name, B_out{ii});

    B_out{ii} = Pipeline_Asmbly_append_preA(B_out{ii});
% end
    % Summary plots
    %                 Pipline_Asmbly_plot(A_out{ii}, [fig_dir filesep method]);
    %     Pipline_Asmbly_plot(P_out{ii}, [fig_dir filesep method filesep 'place']);

    Pipline_Asmbly_plot(B_out{ii}, [fig_dir filesep 'bin_' num2str(bin_size) filesep method filesep 'best']);

    Pipline_Asmbly_plot_SWS(B_out{ii}, [fig_dir filesep 'bin_' num2str(bin_size) filesep method filesep 'best_SWS']);

    
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
    
    if anx_idx(ii)== 0 && novel_idx(ii)==1
        %         B_out{ii} = Pipeline_Asmbly_append_preA(B_out{ii});
        
    end
    
end
novel_idx = logical(novel_idx);
anx_idx = logical(anx_idx);
HS_idx = logical(HS_idx);

lt1_idx = novel_idx & ~anx_idx & ~HS_idx;
lt5_idx = ~novel_idx & ~anx_idx & ~HS_idx;

H1_idx = novel_idx & anx_idx & ~HS_idx;
H5_idx = ~novel_idx & anx_idx & ~HS_idx;

A_all = A_out;
A_out = A_out(1:11);
% % if ~isempty(A_out)
% %     save(['C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\inter\A_out_' method '.mat'], 'A_out')
% % end
% % if ~isempty(P_out)
% save([main_dir  strrep('Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\inter\B_out_', '\', filesep) method '.mat'], 'B_out')
% % end

%% load data that has been processed.

method = 'binary';

load([main_dir  strrep('Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\inter\B_out_', '\', filesep) method '.mat'], 'B_out')

% rename due to progressive structure naming convention above.
A_out = B_out;
clear B_out

% remove pv1254 due to inconsistent running.


exclude_mouse = {'pv1254'};

for ii = length(A_out):-1:1
    
    if contains(A_out{ii}{1}.info.subject, exclude_mouse)
        fprintf('Removing sesson: <strong>%s</strong>\n', A_out{ii}{1}.info.subject);
        rm_idx(ii) = true;
    else
        rm_idx(ii) = false;
    end
    
end

A_out(rm_idx) = [];


novel_idx = []; anx_idx = []; HS_idx = [];
for iA = size(A_out,2):-1:1
    this_f = A_out{iA}{1}.info.session;
    
    %A_out{iA} = Pipeline_Asmbly_append_preA(A_out{iA});
    %A_out{iA} = Pipeline_Asmbly_append_postA(A_out{iA});
    
    if contains(this_f, 'HATDS')
        HS_idx(iA) = 1;
    else
        HS_idx(iA) = 0;
    end
    
    if contains(this_f, 'D1') %|| ~isempty(strfind(f_list(ii).name, 'HATDS'))
        novel_idx(iA) = 1;
    end
    
    if contains(this_f, 'D5')
        novel_idx(iA) = 0;
    end
    
    if contains(this_f, 'HAT')
        anx_idx(iA) = 1;
    else
        anx_idx(iA) = 0;
    end
    
end



%% collect the data

Pre_n_Asmbly = []; Post_n_Asmbly = [];
Pre_r_Asmbly = []; Post_r_Asmbly = [];

S_Pre_n_Asmbly = []; S_Post_n_Asmbly = [];
S_Pre_r_Asmbly = []; S_Post_r_Asmbly = [];

for iB = length(bin_size):-1:1
    
    for iA = size(A_out,2):-1:1
        
        wake_n_Asmbly(iA, iB) = size(A_out{iA}{iB}.P_proj,1);
        
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

s_test = 'ttest';

max_n_A = max([(Pre_n_Asmbly) ; (Post_n_Asmbly) ]);
max_n_A = max_n_A*1.5;
max_r_A = max([(Pre_r_Asmbly) ; (Post_r_Asmbly) ]);
max_r_A = max_r_A*1.5;

max_n_SA = max([(S_Pre_n_Asmbly) ; (S_Post_n_Asmbly) ]);
max_n_SA = max_n_SA*1.5;
max_r_SA = max([(S_Pre_r_Asmbly) ; (S_Post_r_Asmbly) ]);
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


for iB = length(bin_size):-1:1
    n = 4; m = 5;
    figure(100 + iB)
    
    subplot(n,m,1)
    cla
    MS_bar_w_err(Pre_n_Asmbly(n_idx,iB), Post_n_Asmbly(n_idx,iB),c_ord(1,:), 1, s_test);
    
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title(['Novel (' num2str(bin_size(iB)) 's bins)'])
    ylabel('# assemblies')
    ylim([0 max_n_A(iB)])
    
    
    subplot(n,m,2)
    cla
    MS_bar_w_err(Pre_n_Asmbly(f_idx,iB), Post_n_Asmbly(f_idx,iB),c_ord(2,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Fam')
    ylim([0 max_n_A(iB)])
    
    
    subplot(n,m,3)
    cla
    MS_bar_w_err(Pre_n_Asmbly(~a_idx,iB), Post_n_Asmbly(~a_idx,iB),c_ord(3,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Linear')
    ylim([0 max_n_A(iB)])
    
    subplot(n,m,4)
    cla
    MS_bar_w_err(Pre_n_Asmbly(a_idx,iB), Post_n_Asmbly(a_idx,iB),c_ord(4,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Anxiety')
    ylim([0 max_n_A(iB)])
    
    subplot(n,m,5)
    cla
    MS_bar_w_err(Pre_n_Asmbly(HS_idx,iB), Post_n_Asmbly(HS_idx,iB),c_ord(5,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('A. Switch')
    ylim([0 max_n_A(iB)])
    
    
    % same for the rate
    
    subplot(n,m,m+1)
    cla
    MS_bar_w_err(Pre_r_Asmbly(n_idx,iB), Post_r_Asmbly(n_idx,iB),c_ord(1,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title(['Novel (' num2str(bin_size(iB)) 's bins)'])
    ylabel('ReAct/min')
    ylim([0 max_r_A(iB)])
    
    
    subplot(n,m,m+2)
    cla
    MS_bar_w_err(Pre_r_Asmbly(f_idx,iB), Post_r_Asmbly(f_idx,iB),c_ord(2,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Fam')
    ylim([0 max_r_A(iB)])
    
    
    subplot(n,m,m+3)
    cla
    MS_bar_w_err(Pre_r_Asmbly(~a_idx,iB), Post_r_Asmbly(~a_idx,iB),c_ord(3,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Linear')
    ylim([0 max_r_A(iB)])
    
    
    subplot(n,m,m+4)
    cla
    MS_bar_w_err(Pre_r_Asmbly(a_idx,iB), Post_r_Asmbly(a_idx,iB),c_ord(4,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Anxiety')
    ylim([0 max_r_A(iB)])
    
    
    subplot(n,m,m+5)
    cla
    MS_bar_w_err(Pre_r_Asmbly(HS_idx,iB), Post_r_Asmbly(HS_idx,iB),c_ord(5,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('A. Switch')
    ylim([0 max_r_A(iB)])
    
    
    
    % repeat for assemblies with coherent spatial tuning.
    subplot(n,m,(m*2)+1)
    cla
    MS_bar_w_err(S_Pre_n_Asmbly(n_idx,iB), S_Post_n_Asmbly(n_idx,iB),c_ord_s(1,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title(['Novel (' num2str(bin_size(iB)) 's bins)'])
    ylabel({'Spatial A only'; '# assemblies'})
    ylim([0 max_n_SA(iB)])
    
    
    subplot(n,m,(m*2)+2)
    cla
    MS_bar_w_err(S_Pre_n_Asmbly(f_idx,iB), S_Post_n_Asmbly(f_idx,iB),c_ord_s(2,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Fam')
    ylim([0 max_n_SA(iB)])
    
    
    subplot(n,m,(m*2)+3)
    cla
    MS_bar_w_err(S_Pre_n_Asmbly(~a_idx,iB), S_Post_n_Asmbly(~a_idx,iB),c_ord_s(3,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Linear')
    ylim([0 max_n_SA(iB)])
    
    subplot(n,m,(m*2)+4)
    cla
    MS_bar_w_err(S_Pre_n_Asmbly(a_idx,iB), S_Post_n_Asmbly(a_idx,iB),c_ord_s(4,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Anxiety')
    ylim([0 max_n_SA(iB)])
    
    subplot(n,m,(m*2)+5)
    cla
    MS_bar_w_err(S_Pre_n_Asmbly(HS_idx,iB), S_Post_n_Asmbly(HS_idx,iB),c_ord_s(5,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('A. Switch')
    ylim([0 max_n_SA(iB)])
    
    
    % same for the rate
    
    subplot(n,m,(m*3)+1)
    cla
    MS_bar_w_err(S_Pre_r_Asmbly(n_idx,iB), S_Post_r_Asmbly(n_idx,iB),c_ord_s(1,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title(['Novel (' num2str(bin_size(iB)) 's bins)'])
    ylabel({'Spatial A only';'ReAct/min'})
    ylim([0 max_r_SA(iB)])
    
    
    
    subplot(n,m,(m*3)+2)
    cla
    MS_bar_w_err(S_Pre_r_Asmbly(f_idx,iB), S_Post_r_Asmbly(f_idx,iB),c_ord_s(2,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Fam')
    ylim([0 max_r_SA(iB)])
    
    
    subplot(n,m,(m*3)+3)
    cla
    MS_bar_w_err(S_Pre_r_Asmbly(~a_idx,iB), S_Post_r_Asmbly(~a_idx,iB),c_ord_s(3,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Linear')
    ylim([0 max_r_SA(iB)])
    
    
    subplot(n,m,(m*3)+4)
    cla
    MS_bar_w_err(S_Pre_r_Asmbly(a_idx,iB), S_Post_r_Asmbly(a_idx,iB),c_ord_s(4,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('Anxiety')
    ylim([0 max_r_SA(iB)])
    
    
    subplot(n,m,(m*3)+5)
    cla
    MS_bar_w_err(S_Pre_r_Asmbly(HS_idx,iB), S_Post_r_Asmbly(HS_idx,iB),c_ord_s(5,:), 1, s_test);
    set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
    title('A. Switch')
    ylim([0 max_r_SA(iB)])
    
    
end

fprintf('<strong> %s</strong> <strong> %0.02f</strong>\x00B1%0.2f wake assemblies | <strong> %0.02f</strong>\x00B1%0.2f pre | <strong> %0.02f</strong>\x00B1%0.2f post\n',...
    'overall', nanmean(wake_n_Asmbly(:, iB)), std(wake_n_Asmbly(:,iB)), nanmean(Pre_n_Asmbly(:, iB)), std(Pre_n_Asmbly(:,iB)), nanmean(Post_n_Asmbly(:, iB)), std(Post_n_Asmbly(:,iB)))

fprintf('<strong> %s</strong> <strong> %0.02f</strong>\x00B1%0.2f wake assemblies | <strong> %0.02f</strong>\x00B1%0.2f pre | <strong> %0.02f</strong>\x00B1%0.2f post\n',...
    'Novel', nanmean(wake_n_Asmbly(n_idx, iB)), std(wake_n_Asmbly(n_idx,iB)), nanmean(Pre_n_Asmbly(n_idx, iB)), std(Pre_n_Asmbly(n_idx,iB)), nanmean(Post_n_Asmbly(n_idx, iB)), std(Post_n_Asmbly(n_idx,iB)))

fprintf('<strong> %s</strong> <strong> %0.02f</strong>\x00B1%0.2f wake assemblies | <strong> %0.02f</strong>\x00B1%0.2f pre | <strong> %0.02f</strong>\x00B1%0.2f post\n',...
    'Familiar', nanmean(wake_n_Asmbly(~n_idx, iB)), std(wake_n_Asmbly(~n_idx,iB)), nanmean(Pre_n_Asmbly(~n_idx, iB)), std(Pre_n_Asmbly(~n_idx,iB)), nanmean(Post_n_Asmbly(~n_idx, iB)), std(Post_n_Asmbly(~n_idx,iB)))

fprintf('<strong> %s</strong> <strong> %0.02f</strong>\x00B1%0.2f wake assemblies | <strong> %0.02f</strong>\x00B1%0.2f pre | <strong> %0.02f</strong>\x00B1%0.2f post\n',...
    'LT', nanmean(wake_n_Asmbly(~a_idx, iB)), std(wake_n_Asmbly(~a_idx,iB)), nanmean(Pre_n_Asmbly(~a_idx, iB)), std(Pre_n_Asmbly(~a_idx,iB)), nanmean(Post_n_Asmbly(~a_idx, iB)), std(Post_n_Asmbly(~a_idx,iB)))

fprintf('<strong> %s</strong> <strong> %0.02f</strong>\x00B1%0.2f wake assemblies | <strong> %0.02f</strong>\x00B1%0.2f pre | <strong> %0.02f</strong>\x00B1%0.2f post\n',...
    'HAT', nanmean(wake_n_Asmbly(a_idx, iB)), std(wake_n_Asmbly(a_idx,iB)), nanmean(Pre_n_Asmbly(a_idx, iB)), std(Pre_n_Asmbly(a_idx,iB)), nanmean(Post_n_Asmbly(a_idx, iB)), std(Post_n_Asmbly(a_idx,iB)))


%% collect  control reactivation strength

A_Pre_n = []; A_Post_n = [];
A_Pre_r = []; A_Post_r = [];
A_cent = []; A_peak = [];
A_ReAct_all = [];

SA_Pre_n = []; SA_Post_n = [];
SA_Pre_r = []; SA_Post_r = [];
SA_cent = []; SA_peak = [];

A_ReAct = []; J_ReAct = [];
SA_ReAct = []; SJ_ReAct = [];

A_sub_list = []; J_sub_list = [];


for iB = length(bin_size):-1:1
    %     all_cent =
    for iA = size(A_out,2):-1:1
        
        % wake
        
        
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


%% plot  novel VS familiar
% A_out = A_all(1:11);

HS_idx = logical(HS_idx(1:length(A_out)));
n_idx = logical(novel_idx(1:length(A_out)));
a_idx = logical(anx_idx(1:length(A_out)));
lt1_idx = logical(lt1_idx(1:length(A_out)));
lt5_idx = logical(lt5_idx(1:length(A_out)));
H1_idx = logical(H1_idx(1:length(A_out)));
H5_idx = logical(H5_idx(1:length(A_out)));


n_idx = n_idx & ~HS_idx;
a_idx = a_idx & ~HS_idx;
f_idx = ~n_idx & ~HS_idx;


% wider colours
p_cord = MS_linspecer(16);
c_ord = MS_linspecer(5);

figure(4000)
clf
y_lim = [-.12 .12];

subplot(2,4,1)
data_a = A_ReAct(lt1_idx | H1_idx, iB);
data_b = A_ReAct(lt5_idx | H5_idx, iB);

hb = bar([nanmean(data_a), nanmean(data_b)]', 'FaceColor', c_ord(1,:), 'EdgeColor', 'k');
hb.FaceColor = 'flat';
hb.CData(2,:) = c_ord(2,:);
hold on
eb = errorbar([nanmean(data_a), nanmean(data_b)], [MS_SEM(data_a) ,MS_SEM(data_b)]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'novel', 'familiar'}, 'XTickLabelRotation', 45)
ylabel('ReActivation Strength')
ylim(y_lim);
title('Control')

subplot(2,4,2)
data_a = A_ReAct(lt1_idx | lt5_idx, iB);
data_b = A_ReAct(H1_idx | H5_idx, iB);
data_c = A_ReAct(HS_idx, iB);

hb = bar([nanmean(data_a), nanmean(data_b),nanmean(data_c)]', 'FaceColor', c_ord(1,:), 'EdgeColor', 'k');
hb.FaceColor = 'flat';
hb.CData(2,:) = c_ord(2,:);
hb.CData(3,:) = c_ord(4,:);
hold on
eb = errorbar([nanmean(data_a), nanmean(data_b), nanmean(data_c)], [MS_SEM(data_a) ,MS_SEM(data_b), MS_SEM(data_c)]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:3, 'XTickLabel', {'Linear', 'anxiety', 'switch'}, 'XTickLabelRotation', 45)
ylim(y_lim);
title('Control')

% across conditions
subplot(2,4,3)
data_a = A_ReAct(lt1_idx, iB);
data_b = A_ReAct(lt5_idx, iB);
data_c = A_ReAct(H1_idx, iB);
data_d = A_ReAct(H5_idx, iB);
data_e = A_ReAct(HS_idx, iB);


hb = bar([nanmean(data_a), nanmean(data_b),nanmean(data_c) nanmean(data_d),nanmean(data_e)]', 'FaceColor', p_cord(1,:), 'EdgeColor', 'k');
hb.FaceColor = 'flat';
hb.CData(2,:) = p_cord(2,:);
hb.CData(3,:) = p_cord(end-2,:);
hb.CData(4,:) = p_cord(end-1,:);
hb.CData(5,:) = p_cord(end-5,:);
hold on
eb = errorbar([nanmean(data_a), nanmean(data_b), nanmean(data_c) nanmean(data_d),nanmean(data_e)], [MS_SEM(data_a) ,MS_SEM(data_b), MS_SEM(data_c), MS_SEM(data_d), MS_SEM(data_e)]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:5, 'XTickLabel', {'LT1', 'LT5', 'H1' , 'H5', 'switch'}, 'XTickLabelRotation', 45)
ylim(y_lim);
title('Control')



% Square_subplots

SetFigure([],gcf)

saveas(gcf,[fig_dir filesep 'ReAct_summary_p' num2str(bin_size(iB)) '.png'])

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
        A_n_sub = [A_n_sub repmat(ii, 1, length(A_ReAct_all{iB}{ii}))];
    end
    
    if ismember(ii, find(f_idx))
        A_f_cent = [A_f_cent A_cent{iB}{ii}];
        A_f_peak = [A_f_peak A_peak{iB}{ii}];
        A_f_ReAct = [A_f_ReAct A_ReAct_all{iB}{ii}'];
        A_f_sub = [A_f_sub repmat(ii, 1, length(A_ReAct_all{iB}{ii}))];
    end
    
end

% control
figure(5000)
subplot(2,2,1)
scatter(A_n_cent, A_n_ReAct, 35, A_n_sub, 'filled');
[c, p] = corrcoef(A_n_cent, A_n_ReAct);
text(-1, .7, ['r = ' num2str(c(1,2),2) ' p= ' num2str(p(1,2),2)], 'FontWeight', 'bold')
xlabel('(<- coherent) Z Map centroid')
ylabel('ReAct Str.')
title('Control novel (Color by subject)')
ylim([-1 1]);
xlim([-5 5]);
lsline

subplot(2,2,2)
scatter(A_f_cent, A_f_ReAct, 35, A_f_sub, 'filled', 'd');
[c, p] = corrcoef(A_f_cent, A_f_ReAct);
text(-1, .7, ['r = ' num2str(c(1,2),2) ' p= ' num2str(p(1,2),2)], 'FontWeight', 'bold')
xlabel('(<- coherent) Z Map centroid')
ylabel('ReAct Str.')
title('familiar')
ylim([-1 1]);
xlim([-5 5]);
lsline

subplot(2,2,3)
scatter(A_n_peak, A_n_ReAct, 35, A_n_sub, 'filled');
[c, p] = corrcoef(A_n_peak, A_n_ReAct);
text(-1, .7, ['r = ' num2str(c(1,2),2) ' p= ' num2str(p(1,2),2)], 'FontWeight', 'bold')
ylabel('ReAct Str.')
xlabel('(<- coherent) Z Map peak')
ylim([-1 1]);
xlim([-5 5]);
lsline

subplot(2,2,4)
scatter(A_f_peak, A_f_ReAct, 35, A_f_sub, 'filled', 'd');
[c, p] = corrcoef(A_f_peak, A_f_ReAct);
text(-1, .7, ['r = ' num2str(c(1,2),2) ' p= ' num2str(p(1,2),2)], 'FontWeight', 'bold')
xlabel('(<- coherent) Z Map peak')
ylabel('ReAct Str.')
ylim([-1 1]);
xlim([-5 5]);
lsline


%% Collect the change in significant CCF
pre_sig_cff = []; post_sig_cff = [];
for iB = length(bin_size):-1:1
    
    for iS = size(A_out,2):-1:1
        
        pre_sig_cff(iS, iB) = A_out{iS}{iB}.REM_Pre_psig_cff;
        post_sig_cff(iS, iB) = A_out{iS}{iB}.REM_Post_psig_cff;
    end
    
    
    figure(10+iB)
    %     c_ord = MS_linspecer(4);
    ax(1) = subplot(2,5,1);
    cla
    MS_bar_w_err(pre_sig_cff(n_idx,iB), post_sig_cff(n_idx,iB), c_ord(1,:), 1, s_test);
    set(gca, 'XTickLabel', {'pre', 'post'})
    title('novel');
    %     title(['CCF sig (bin: ' num2str(bin_size(iB)) ')'])
    ylabel({'Cross correlation'; '% significant of wake'})
    
    ax(2) = subplot(2,5,2);
    cla
    MS_bar_w_err(pre_sig_cff(f_idx,iB), post_sig_cff(f_idx,iB), c_ord(2,:), 1, s_test);
    set(gca, 'XTickLabel', {'pre', 'post'})
    title('familiar');
    
    
    ax(3) = subplot(2,5,3);
    cla
    MS_bar_w_err(pre_sig_cff(lt1_idx | lt5_idx,iB), post_sig_cff(lt1_idx | lt5_idx,iB), c_ord(3,:), 1, s_test);
    set(gca, 'XTickLabel', {'pre', 'post'})
    title('Linear');
    
    
    ax(4) = subplot(2,5,4);
    cla
    MS_bar_w_err(pre_sig_cff(H1_idx | H5_idx,iB), post_sig_cff(H1_idx | H5_idx,iB), c_ord(4,:), 1, s_test);
    set(gca, 'XTickLabel', {'pre', 'post'})
    title('Anxiety');
    
    ax(5) =  subplot(2,5,5);
    cla
    MS_bar_w_err(pre_sig_cff(H1_idx | H5_idx,iB), post_sig_cff(H1_idx | H5_idx,iB), c_ord(5,:), 1, s_test);
    set(gca, 'XTickLabel', {'pre', 'post'})
    title('switch');
    
    linkprop(ax, 'ylim')
    ylim([0 max([pre_sig_cff,post_sig_cff], [], 'all')*1.1])
    %     Square_subplots;
    SetFigure([], gcf)
    %     title('J20')
    saveas(gcf,[fig_dir filesep 'cff_summary_p' num2str(strrep(bin_size(iB), '.', 'p')) '.png'])
    
end


%% qualitfy reactivation saptial biases

data_in = A_out;
Pre_hist = []; Post_hist = [];
A_hist_pre = []; A_hist_post = [];

for kk = 1
    
    data_in = A_out;
    
    
    for iB = length(data_in{ii}):-1:1
        
        all_cnt_prct = []; Pre_cnt_prct = []; Post_cnt_prct = [];
        for ii = length(data_in):-1:1
            
            % data shorthand
            this_A = data_in{ii}{iB};
            
            % make arrays to accrue reactivation locations;
            Pre_ReAct_mat = []; Post_ReAct_mat = []; Wake_ReAct_mat = [];
            
            cent_z = []; peak_z = []; M_cent = [];
            
            Pre_cnt = sum((this_A.REM_Pre_proj > this_A.REM_Pre_stats.R_thresh), 2);
            Post_cnt = sum((this_A.REM_Post_proj > this_A.REM_Post_stats.R_thresh),2);
            
            Wake_cnt = sum((this_A.P_proj > 10), 2);
            
            
            
            for iA = length(this_A.P_pos):-1:1
                
                cent_z(iA) = this_A.map{iA}.cent_z;
                peak_z(iA) = this_A.map{iA}.peak_z;
                
                M_cent(iA) = mean(this_A.map{iA}.cent);
                if cent_z(iA) < -1.96
                    Pre_ReAct_mat = [Pre_ReAct_mat repmat(M_cent(iA),1, Pre_cnt(iA))];
                    Post_ReAct_mat = [Post_ReAct_mat repmat(M_cent(iA),1, Post_cnt(iA))];
                    Wake_ReAct_mat = [Wake_ReAct_mat repmat(M_cent(iA),1, Wake_cnt(iA))];
                    
                end
            end
            all_cnt_prct(ii) = sum(cent_z < -1.96) / length(cent_z);
            Pre_cnt_prct(ii) = sum((cent_z < -1.96) & (this_A.REM_Pre_stats.p_val < 0.05))/sum(cent_z < -1.96);
            Post_cnt_prct(ii) = sum((cent_z < -1.96) & (this_A.REM_Post_stats.p_val < 0.05))/sum(cent_z < -1.96);
            disp(max(Pre_ReAct_mat))
            
            p_bins = 0:5:100;
            [Pre_hist{ii, iB}] = histcounts(Pre_ReAct_mat*2.5, p_bins);
            [Post_hist{ii, iB}] = histcounts(Post_ReAct_mat*2.5, p_bins);
            [Wake_hist{ii, iB}] = histcounts(Wake_ReAct_mat*2.5, p_bins);
            
            if ~isempty(Pre_hist{ii, iB})
                Pre_hist{ii, iB} = Pre_hist{ii, iB}./((this_A.REM_Pre_tvec(end) - this_A.REM_Pre_tvec(1))/60);
            else
                Pre_hist{ii, iB} = NaN(1,length(p_bins(1:end)));
            end
            
            if ~isempty(Post_hist{ii, iB})
                Post_hist{ii, iB} = Post_hist{ii, iB}./((this_A.REM_Post_tvec(end) - this_A.REM_Post_tvec(1))/60);
            else
                Post_hist{ii, iB} = NaN(1,length(p_bins(1:end)));
            end
            
            if kk == 1
                A_hist_pre(ii,:) = Pre_hist{ii,iB};
                A_hist_post(ii,:) = Post_hist{ii,iB};
                A_hist_wake(ii,:) = Wake_hist{ii, iB};
            end
            
        end
        
    end
    
    if kk == 1
        A_cnt_prct = all_cnt_prct;
        A_Pre_cnt_prct = Pre_cnt_prct;
        A_Post_cnt_prct = Post_cnt_prct;
    end
end

p_centr = p_bins(1:end-1)+3/2;
a_ord = MS_linspecer(2);

figure(300)
clf
subplot(5,1,1)
hold on
bar(p_centr, nanmean(A_hist_pre(lt1_idx,:))./max(nanmean(A_hist_pre(lt1_idx,:))), 1, 'facecolor', a_ord(2,:), 'FaceAlpha', .3);
bar(p_centr, nanmean(A_hist_post(lt1_idx,:))./max(nanmean(A_hist_post(lt1_idx,:))),1, 'facecolor', a_ord(1,:), 'FaceAlpha', .3);
set(gca, 'XTickLabel', [], 'YTick', [0 1]);
legend({'pre', 'post'}, 'box', 'off', 'Orientation', 'horizontal', 'Location', 'northwest')


subplot(5,1,2)
hold on
bar(p_centr, nanmean(A_hist_pre(lt5_idx,:))./max(nanmean(A_hist_pre(lt5_idx,:))), 'facecolor', a_ord(2,:), 'FaceAlpha', .3);
bar(p_centr, nanmean(A_hist_post(lt5_idx,:))./max(nanmean(A_hist_post(lt5_idx,:))), 'facecolor', a_ord(1,:), 'FaceAlpha', .3);
set(gca, 'XTickLabel', [], 'YTick', [0 1]);

subplot(5,1,3)
hold on
bar(p_centr, nanmean(A_hist_pre(H1_idx,:))./max(nanmean(A_hist_pre(H1_idx,:))), 'facecolor', a_ord(2,:), 'FaceAlpha', .3);
bar(p_centr, nanmean(A_hist_post(H1_idx,:))./max(nanmean(A_hist_post(H1_idx,:))), 'facecolor', a_ord(1,:), 'FaceAlpha', .3);
set(gca, 'XTickLabel', [], 'YTick', [0 1]);


subplot(5,1,4)
hold on
bar(p_centr, nanmean(A_hist_pre(H5_idx,:))./max(nanmean(A_hist_pre(H5_idx,:))), 'facecolor', a_ord(2,:), 'FaceAlpha', .3);
bar(p_centr, nanmean(A_hist_post(H5_idx,:))./max(nanmean(A_hist_post(H5_idx,:))), 'facecolor', a_ord(1,:), 'FaceAlpha', .3);
set(gca, 'XTickLabel', [], 'YTick', [0 1]);

subplot(5,1,5)
hold on
bar(p_centr, nanmean(A_hist_pre(HS_idx,:))./max(nanmean(A_hist_pre(HS_idx,:))), 'facecolor', a_ord(2,:), 'FaceAlpha', .3);
bar(p_centr, nanmean(A_hist_post(HS_idx,:))./max(nanmean(A_hist_post(HS_idx,:))), 'facecolor', a_ord(1,:), 'FaceAlpha', .3);
% stem(p_centr, nanmean(A_hist_pre(HS_idx,:))./max(nanmean(A_hist_pre(HS_idx,:))), 'color', a_ord(2,:));
% stem(p_centr, nanmean(A_hist_post(HS_idx,:))./max(nanmean(A_hist_post(HS_idx,:))), 'color', a_ord(1,:));
set(gca, 'YTick', [0 1], 'XTick', [0 100]);
xlabel('position on track (cm)')



% figure(301)
% clf
% subplot(3,1,1)
% hold on
% area(p_centr, nanmean(J20_hist_pre(J20_novel_idx,:)), 'facecolor', c_ord(3,:), 'FaceAlpha', .6);
% area(p_centr, nanmean(J20_hist_post(J20_novel_idx,:)), 'facecolor', c_ord(4,:), 'FaceAlpha', .6);
%
% subplot(3,1,2)
% hold on
% area(p_centr, nanmean(J20_hist_pre(D3_idx,:)), 'facecolor', c_ord(3,:), 'FaceAlpha', .6);
% area(p_centr, nanmean(J20_hist_post(D3_idx,:)), 'facecolor', c_ord(4,:), 'FaceAlpha', .6);
% subplot(3,1,3)
% hold on
% area(p_centr, nanmean(J20_hist_pre(D5_idx,:)), 'facecolor', c_ord(3,:), 'FaceAlpha', .6);
% area(p_centr, nanmean(J20_hist_post(D5_idx,:)), 'facecolor', c_ord(4,:), 'FaceAlpha', .6);
%
%
% % bar(nanmean([Post_hist{1,1}; Post_hist{4,1}; Post_hist{7,1}])


%% quantify using whole spatial map averaging

data_in = A_out;
Diff_ReAct_map = [];
Pre_ReAct_mean_map = [];
Post_ReAct_mean_map = [];

% loop over window size
for iB = length(data_in{ii}):-1:1
    
    all_cnt_prct = []; Pre_cnt_prct = []; Post_cnt_prct = [];
    % loop over sessions
    for ii = length(data_in):-1:1
        
        % data shorthand
        this_A = data_in{ii}{iB};
        
        % make arrays to accrue reactivation locations;
        Pre_ReAct_mat = []; Post_ReAct_mat = []; Wake_ReAct_mat = [];
        Pre_ReAct_map = cell(length(data_in),1); Post_ReAct_map = cell(length(data_in),1); Wake_ReAct_map = cell(length(data_in),1);
        
        cent_z = []; peak_z = []; M_cent = [];
        
        Pre_cnt = sum((this_A.REM_Pre_proj > this_A.REM_Pre_stats.R_thresh), 2);
        Post_cnt = sum((this_A.REM_Post_proj > this_A.REM_Post_stats.R_thresh),2);
        
        Wake_cnt = sum((this_A.P_proj > 10), 2);
        
        
        
        for iA = length(this_A.P_pos):-1:1
            
            cent_z(iA) = this_A.map{iA}.cent_z;
            peak_z(iA) = this_A.map{iA}.peak_z;
            
            M_cent(iA) = mean(this_A.map{iA}.cent);
            
            if cent_z(iA) < -1.96
                Pre_ReAct_map{ii} = [Pre_ReAct_map{ii} ;repmat(this_A.map{iA}.map_mean./max(this_A.map{iA}.map_mean), Pre_cnt(iA),1)];
                Post_ReAct_map{ii} = [Post_ReAct_map{ii} ;repmat(this_A.map{iA}.map_mean./max(this_A.map{iA}.map_mean), Post_cnt(iA),1)];
                Wake_ReAct_map{ii} = [Wake_ReAct_map{ii} ;repmat(this_A.map{iA}.map_mean./max(this_A.map{iA}.map_mean), Wake_cnt(iA),1)];
                
            end
        end
        
        if isempty(Pre_ReAct_map{ii})
            Pre_ReAct_map{ii} = NaN(1,length(p_bins(1:end)));
        end
        
        if isempty(Post_ReAct_map{ii})
            Post_ReAct_map{ii} = NaN(1,length(p_bins(1:end)));
        end
        
        if isempty(Wake_ReAct_map{ii})
            Wake_ReAct_map{ii} = NaN(1,length(p_bins(1:end)));
        end
        
        Pre_ReAct_mean_map(ii,:) = mean(Pre_ReAct_map{ii});
        Post_ReAct_mean_map(ii,:) = mean(Post_ReAct_map{ii});
        Diff_ReAct_map(ii,:) = Pre_ReAct_mean_map(ii,:)./Post_ReAct_mean_map(ii,:);
        Wake_ReAct_mean_map(ii,:) = mean(Wake_ReAct_map{ii});
        
        
    end
    
end
%
a_ord = MS_linspecer(2);

figure(300)
clf
subplot(5,1,1)
hold on
bar(this_A.map{1}.bins, nanmean(Pre_ReAct_mean_map(lt1_idx,:)), 1, 'facecolor', a_ord(2,:), 'FaceAlpha', .3);
bar(this_A.map{1}.bins, nanmean(Post_ReAct_mean_map(lt1_idx,:)),1, 'facecolor', a_ord(1,:), 'FaceAlpha', .3);
plot(this_A.map{1}.bins, nanmean(Wake_ReAct_mean_map(lt1_idx,:)), 'color', 'k');
legend({'pre', 'post', 'wake'}, 'box', 'off', 'Orientation', 'horizontal', 'Location', 'northwest')

set(gca, 'XTickLabel', [], 'YTick', [0 1]);

subplot(5,1,2)
hold on
bar(this_A.map{1}.bins, nanmean(Pre_ReAct_mean_map(lt5_idx,:)), 1, 'facecolor', a_ord(2,:), 'FaceAlpha', .3);
bar(this_A.map{1}.bins, nanmean(Post_ReAct_mean_map(lt5_idx,:)),1, 'facecolor', a_ord(1,:), 'FaceAlpha', .3);
plot(this_A.map{1}.bins, nanmean(Wake_ReAct_mean_map(lt5_idx,:)), 'color', 'k');

set(gca, 'XTickLabel', [], 'YTick', [0 1]);

subplot(5,1,3)
hold on
rectangle('position', [(this_A.map{1}.bins(end) - this_A.map{1}.bins(1) +this_A.map{1}.bin_size)/2, 0, (this_A.map{1}.bins(end) - this_A.map{1}.bins(1) +this_A.map{1}.bin_size)/2, max([nanmean(Pre_ReAct_mean_map(H1_idx,:)) nanmean(Post_ReAct_mean_map(H1_idx,:))])*1.3], ...
    'FaceColor', [.8 .8 .8], 'EdgeColor', [1 1 1])
bar(this_A.map{1}.bins, nanmean(Pre_ReAct_mean_map(H1_idx,:)), 1, 'facecolor', a_ord(2,:), 'FaceAlpha', .3);
bar(this_A.map{1}.bins, nanmean(Post_ReAct_mean_map(H1_idx,:)),1, 'facecolor', a_ord(1,:), 'FaceAlpha', .3);
plot(this_A.map{1}.bins, nanmean(Wake_ReAct_mean_map(H1_idx,:)), 'color', 'k');
xlim([this_A.map{1}.bins(1)-.5 this_A.map{1}.bins(end)+0.5])
ylim([0 max([nanmean(Pre_ReAct_mean_map(H1_idx,:)) nanmean(Post_ReAct_mean_map(H1_idx,:))])*1.1])
set(gca, 'XTickLabel', [], 'YTick', [0 1]);


subplot(5,1,4); cla
hold on
rectangle('position', [(this_A.map{1}.bins(end) - this_A.map{1}.bins(1) +this_A.map{1}.bin_size)/2, 0, (this_A.map{1}.bins(end) - this_A.map{1}.bins(1) +this_A.map{1}.bin_size)/2, max([nanmean(Pre_ReAct_mean_map(H5_idx,:)) nanmean(Post_ReAct_mean_map(H5_idx,:))])*1.3], ...
    'FaceColor', [.8 .8 .8], 'EdgeColor', [1 1 1])
bar(this_A.map{1}.bins, nanmean(Pre_ReAct_mean_map(H5_idx,:)), 1, 'facecolor', a_ord(2,:), 'FaceAlpha', .3);
bar(this_A.map{1}.bins, nanmean(Post_ReAct_mean_map(H5_idx,:)),1, 'facecolor', a_ord(1,:), 'FaceAlpha', .3);
plot(this_A.map{1}.bins, nanmean(Wake_ReAct_mean_map(H5_idx,:)), 'color', 'k');
xlim([this_A.map{1}.bins(1)-.5 this_A.map{1}.bins(end)+0.5])
ylim([0 max([nanmean(Pre_ReAct_mean_map(H5_idx,:)) nanmean(Post_ReAct_mean_map(H5_idx,:))])*1.1])
set(gca, 'XTickLabel', [], 'YTick', [0 1]);

subplot(5,1,5)
hold on
rectangle('position', [this_A.map{1}.bins(1) - this_A.map{1}.bin_size/2, 0, (this_A.map{1}.bins(end) - this_A.map{1}.bins(1) +this_A.map{1}.bin_size)/2, max([nanmean(Pre_ReAct_mean_map(HS_idx,:)) nanmean(Post_ReAct_mean_map(HS_idx,:))])*1.3], ...
    'FaceColor', [.8 .8 .8], 'EdgeColor', [1 1 1]);
bar(this_A.map{1}.bins, nanmean(Pre_ReAct_mean_map(HS_idx,:)), 1, 'facecolor', a_ord(2,:), 'FaceAlpha', .3);
bar(this_A.map{1}.bins, nanmean(Post_ReAct_mean_map(HS_idx,:)),1, 'facecolor', a_ord(1,:), 'FaceAlpha', .3);
plot(this_A.map{1}.bins, nanmean(Wake_ReAct_mean_map(HS_idx,:)), 'color', 'k');

xlim([this_A.map{1}.bins(1)-.5 this_A.map{1}.bins(end)+0.5])
ylim([0 max([nanmean(Pre_ReAct_mean_map(HS_idx,:)) nanmean(Post_ReAct_mean_map(HS_idx,:))])*1.1])
set(gca, 'XTickLabel', [], 'YTick', [0 1]);
set(gca, 'YTick', [0 1], 'XTick', [0 100]);
xlabel('position on track (cm)')

%%
figure(10101)
clf
off_set = 0:2.5:10;
p_centr = this_data.map{1}.bins;
for ii = 1:5
    if ii == 1
        this_idx = lt1_idx;
    elseif ii == 2
        this_idx = lt5_idx;
    elseif ii == 3
        this_idx = H1_idx;
    elseif ii == 4
        this_idx = H5_idx;
    elseif ii == 5
        this_idx = HS_idx;
    end
    % add 'walls'
    %      if ii == 1 || ii == 2
    %         rectangle('position', [off_set(ii), p_centr(1), mode(diff(off_set))*.75, p_centr(end)- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    %     elseif ii == 3 || ii ==4
    %         rectangle('position', [off_set(ii), p_centr(ceil(length(p_centr)/2)), mode(diff(off_set))*.75, p_centr(end)- p_centr(ceil(length(p_centr)/2))], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    %     elseif ii ==5
    %         rectangle('position', [off_set(ii), p_centr(1), mode(diff(off_set))*.75, p_centr(ceil(length(p_centr)/2))- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    %        end
    %
    hold on
    norm_val= max([nanmean(Diff_ReAct_map(this_idx,:))]);
    plot(((nanmean(Diff_ReAct_map(this_idx,:))./norm_val)-1)*100+off_set(ii)*100, p_centr,   'linewidth', 2, 'Color', f_ord(ii,:))
    % stem(p_centr,(((nanmean(Diff_ReAct_map(this_idx,:))./norm_val)-1)*100)+off_set(ii)*100,  'linewidth', 2, 'Color', f_ord(ii,:))
    
    xline(off_set(ii)*100)
    if ii == 1
        %         text(off_set(ii)-.1, 97, '\leftarrow', 'HorizontalAlignment','right', 'Interpreter', 'TeX', 'fontsize', 8)
        text((off_set(ii)-.1)*100, p_centr(end), 'pre', 'HorizontalAlignment','right', 'Interpreter', 'TeX', 'fontsize', 8, 'VerticalAlignment', 'baseline')
        %         text(off_set(ii)+.1, 97, '\rightarrow', 'HorizontalAlignment','left', 'Interpreter', 'TeX', 'fontsize', 8)
        text((off_set(ii)+.1)*100, p_centr(end), 'post', 'HorizontalAlignment','left', 'Interpreter', 'TeX', 'fontsize', 8, 'VerticalAlignment', 'baseline')
        
    end
end
% set(gca, 'xtick', off_set, 'ytick', [0 100])
% set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
% ylabel('position (cm)')
% axis('square')
% title('REM assembly tuning')
%% convert to HD5


for ii = length(A_out):-1:1
    
    
    
    fname = ['assembly_' A_out{ii}{1}.info.session '_' A_out{ii}{1}.info.subject '.h5'];
    
    if exist(fname, 'file')
        delete(fname)
    end
    
    hdf5write(fname, '/mouse', string(A_out{ii}{1}.info.subject));
    hdf5write(fname, '/condition', string(A_out{ii}{1}.info.session),'WriteMode', 'append');
    
    this_pre_rem_t = 0:1/33:(length(A_out{ii}{1}.REM_pre_in)/33);
    this_pre_rem_t = this_pre_rem_t(1:end-1);
    
    this_post_rem_t = 0:1/33:(length(A_out{ii}{1}.REM_post_in)/33);
    this_post_rem_t = this_post_rem_t(1:end-1);
    
    pre_React = []; pre_ID = []; pre_sig = []; pre_str = [];
    
    post_React = []; post_ID = []; post_sig = []; post_str = [];
    
    map_loc = []; wake_rate = []; Asmbly_dir = [];
    
    nWake_R = [];
    
    weights = A_out{ii}{1}.P_temp;
    
    A_cell_id = [];
    
    for aa = 1:size(A_out{ii}{1}.REM_Pre_proj,1)
        this_R_t  =  A_out{ii}{1}.REM_Pre_tvec(find(A_out{ii}{1}.REM_Pre_proj(aa, :) > A_out{ii}{1}.REM_Pre_stats.R_thresh));
        pre_React  = [pre_React, (nearest_idx3(this_R_t, this_pre_rem_t))'];
        pre_str = [pre_str, A_out{ii}{1}.REM_Pre_proj(aa, find(A_out{ii}{1}.REM_Pre_proj(aa, :) > A_out{ii}{1}.REM_Pre_stats.R_thresh))];
        
        pre_ID = [pre_ID repmat(aa,1,length(this_R_t))];
        if A_out{ii}{1}.REM_Pre_stats.p_val(aa) < 0.05
            pre_sig = [pre_sig, ones(1, length(this_R_t))];
        else
            pre_sig = [pre_sig, zeros(1, length(this_R_t))];
        end
        
        
        this_R_t  =  A_out{ii}{1}.REM_Post_tvec(find(A_out{ii}{1}.REM_Post_proj(aa, :) > A_out{ii}{1}.REM_Post_stats.R_thresh));
        post_React  = [post_React, (nearest_idx3(this_R_t, this_post_rem_t))'];
        post_str = [post_str, A_out{ii}{1}.REM_Post_proj(aa, find(A_out{ii}{1}.REM_Post_proj(aa, :) > A_out{ii}{1}.REM_Post_stats.R_thresh))];
        
        post_ID = [post_ID repmat(aa,1,length(this_R_t))];
        if A_out{ii}{1}.REM_Post_stats.p_val(aa) < 0.05
            post_sig = [post_sig, ones(1, length(this_R_t))];
        else
            post_sig = [post_sig, zeros(1, length(this_R_t))];
        end
        
        % reactivation strength
        ReAct_post(aa) = mean(A_out{ii}{1}.REM_Post_proj(aa,:)) - mean(A_out{ii}{1}.REM_Pre_proj(aa,:));
        
        % get the mean place field for the assembly
        [~, idx] = max(A_out{ii}{1}.map{aa}.map_mean);
        map_loc(aa) = A_out{ii}{1}.map{aa}.bins(idx);
        
        wake_rate(aa) = length(A_out{ii}{1}.P_loc{aa}.loc);
        
        Asmbly_dir(aa) = mean(A_out{ii}{1}.P_loc{aa}.loc_dir);
        
        % weights
        
        A_cell_id(:,aa) = ismember(1:length(A_out{ii}{1}.P_temp),A_out{ii}{1}.P_pos{aa})';
        
        
        nWake_R(aa) = sum(A_out{ii}{1}.P_proj(aa,:) > 10) ./ (A_out{ii}{1}.wake_tvec(end) - A_out{ii}{1}.wake_tvec(1));
        
    end
    
    ReAct_Str_all_post = mean(ReAct_post);
    
    Total_move = nansum(abs(diff(A_out{ii}{1}.behav.position(:,1))));
    
    prct_move = sum(A_out{ii}{1}.move_idx)./length(A_out{ii}{1}.move_idx);
    
    
    REM_Pre_nSigA = sum(A_out{ii}{1}.REM_Pre_stats.p_val <0.05);
    
    REM_Post_nSigA = sum(A_out{ii}{1}.REM_Post_stats.p_val <0.05);
    
    SWS_Pre_nSigA = sum(A_out{ii}{1}.SWS_Pre_stats.p_val <0.05);
    
    SWS_Post_nSigA = sum(A_out{ii}{1}.SWS_Post_stats.p_val <0.05);
    
    pREM_map = [];
    for pp = length(A_out{ii}{1}.pREM_Place_map):-1:1
        pREM_map = [pREM_map; A_out{ii}{1}.pREM_Place_map{pp}.map];
    end
    
    p_wake_react = [];
    
    % wake assemblies projected from Pre REM data.
    %     for aa = 1:size( A_out{ii}{1}.REM_Wake_proj,1)
    %         this_R_t  =  A_out{ii}{1}.wake_tvec(find(A_out{ii}{1}.REM_Wake_proj(aa, :) > A_out{ii}{1}.REM_Pre_stats.R_thresh));
    %         p_wake_react  = [p_wake_react,this_R_t'];
    %
    %         pre_str = [pre_str, A_out{ii}{1}.REM_Pre_proj(aa, find(A_out{ii}{1}.REM_Pre_proj(aa, :) > A_out{ii}{1}.REM_Pre_stats.R_thresh))];
    %         pre_ID = [pre_ID repmat(aa,1,length(this_R_t))];
    %         if A_out{ii}{1}.REM_Pre_stats.p_val(aa) < 0.05
    %             pre_sig = [pre_sig, ones(1, length(this_R_t))];
    %         else
    %             pre_sig = [pre_sig, zeros(1, length(this_R_t))];
    %         end
    %
    
    %         % get the mean place field for the assembly
    %         [~, idx] = max(A_out{ii}{1}.map{aa}.map_mean);
    %         map_loc(aa) = A_out{ii}{1}.map{aa}.bins(idx);
    %
    %         wake_rate(aa) = length(A_out{ii}{1}.P_loc{aa}.loc);
    %
    %         Asmbly_dir(aa) = mean(A_out{ii}{1}.P_loc{aa}.loc_dir);
    %
    %         % weights
    %
    %         A_cell_id(:,aa) = ismember(1:length(A_out{ii}{1}.P_temp),A_out{ii}{1}.P_pos{aa})';
    %
    
    %     end
    
    
    
    hdf5write(fname, '/pre_rem_A_react_idx', int16(pre_React),'WriteMode', 'append');
    hdf5write(fname, '/pre_rem_A_ID', int8(pre_ID),'WriteMode', 'append');
    hdf5write(fname, '/pre_rem_A_sig', int8(pre_sig),'WriteMode', 'append');
    hdf5write(fname, '/pre_rem_A_str', pre_str,'WriteMode', 'append');
    hdf5write(fname, '/pre_rem_A_sig_react', int8(REM_Pre_nSigA),'WriteMode', 'append');
    
    hdf5write(fname, '/post_rem_react_idx',int16(post_React),'WriteMode', 'append');
    hdf5write(fname, '/post_rem_A_ID', int8(post_ID),'WriteMode', 'append');
    hdf5write(fname, '/post_rem_A_sig', int8(post_sig),'WriteMode', 'append');
    hdf5write(fname, '/post_rem_A_str', post_str,'WriteMode', 'append');
    hdf5write(fname, '/post_rem_A_sig_react', int8(REM_Post_nSigA),'WriteMode', 'append');
    
    
    hdf5write(fname, '/ReAct_Str_all_post', (ReAct_Str_all_post),'WriteMode', 'append');
    
    hdf5write(fname, '/REM_pre_rThresh', A_out{ii}{1}.REM_Pre_stats.R_thresh,'WriteMode', 'append');
    hdf5write(fname, '/REM_post_rThresh', A_out{ii}{1}.REM_Post_stats.R_thresh,'WriteMode', 'append');
    
    hdf5write(fname, '/nWake_R', mean(nWake_R),'WriteMode', 'append');
    
    hdf5write(fname, '/Total_move', mean(Total_move),'WriteMode', 'append');
    hdf5write(fname, '/prct_move', mean(prct_move),'WriteMode', 'append');
    
    
    % sws
    hdf5write(fname, '/pre_sws_A_sig_react', int8(SWS_Pre_nSigA),'WriteMode', 'append');
    hdf5write(fname, '/post_sws_A_sig_react', int8(SWS_Post_nSigA),'WriteMode', 'append');
    
    
    hdf5write(fname, '/wake_rate', int8(wake_rate),'WriteMode', 'append');
    hdf5write(fname, '/map_loc', map_loc,'WriteMode', 'append');
    
    hdf5write(fname, '/wake_a_dir', int8(Asmbly_dir),'WriteMode', 'append');
    
    %     if strcmp(A_out{ii}{1}.info.session, 'LTD1')
    Pre_REM_nA = size(A_out{ii}{1}.pREM_temp,2);
    
    if Pre_REM_nA == 0
        Pre_REM_wAct = 0;
    else
        Pre_REM_wAct = sum(A_out{ii}{1}.pREM_Wake_proj > 8, 2);
    end
    
    fprintf('%s: Pre-ReM_nA = %0.0f | Pre_REM_wAct = %0.0f\n', fname(1:end-3), Pre_REM_nA, sum(Pre_REM_wAct >0))
    
    hdf5write(fname, '/pre_rem_nA', int8(Pre_REM_nA),'WriteMode', 'append');
    hdf5write(fname, '/pre_rem_nWake_A', int8(Pre_REM_wAct),'WriteMode', 'append');
    %     end
    
    %weights
    hdf5write(fname, '/weights', weights,'WriteMode', 'append');
    hdf5write(fname, '/A_sig_cells', int8(A_cell_id),'WriteMode', 'append');
    
    
    
    %
    %     hdf5write(fname, '/wake_proj', (A_out{ii}{1}.P_proj),'WriteMode', 'append');
    %     hdf5write(fname, '/pre_rem_proj', (A_out{ii}{1}.REM_Pre_proj),'WriteMode', 'append');
    %     hdf5write(fname, '/post_rem_proj', (A_out{ii}{1}.REM_Post_proj),'WriteMode', 'append');
    %
    %     hdf5write(fname, '/ReAct_str', (A_out{ii}{1}.ReAct),'WriteMode', 'append');
    %
    %     hdf5write(fname, '/pre_rem_Rthresh', (A_out{ii}{1}.REM_Pre_stats.R_thresh),'WriteMode', 'append');
    %     hdf5write(fname, '/post_rem_Rthresh',(A_out{ii}{1}.REM_Post_stats.R_thresh),'WriteMode', 'append');
    
    
end

%% read it back
MS_h5_to_stuct(fname)


%% get some basic stats
h_idx = dir('assembly*');

exclude_mouse = {'pv1254'};

for ii = length(h_idx):-1:1
    
    if contains(h_idx(ii).name, exclude_mouse)
        fprintf('Removing sesson: <strong>%s</strong>\n', h_idx(ii).name);
        h_idx(ii) = [];
    end
    
end


conds  = {'LTD1 ', 'LTD5 ', 'HATD1 ', 'HATD5 ', 'HATDSwitch '};
mice = {'pv1043', 'pv1060', 'pv1069', 'pv1191', 'pv1192', 'pv1252', 'pv1254'};
Cond = []; Sub = []; nWake_A = [];nWake_R =[];  nPre_A = []; nPre_Act = []; nPre_SWS_A = []; nPost_SWS_A = [];
PreA_str = []; PostA_str= [];
nPre_sig = []; nPost_sig = [];
nPre_A_prct = []; ReAct_str = [];
wake_move = []; wake_prct_move = [];
R_thresh = [];

for hh = length(h_idx):-1:1
    disp(h_idx(hh).name)
    this_h5 = MS_h5_to_stuct(h_idx(hh).name);
    
    Cond(hh) = find(contains(conds, this_h5.condition{1,1}(1:end-1)));
    Sub(hh) = find(contains(mice, this_h5.mouse{1,1}(1:end-1)));
    
    nWake_A(hh) = size(this_h5.A_sig_cells,2);
    nWake_R(hh) = mean(this_h5.nWake_R);
    wake_move(hh) = this_h5.Total_move/100;
    wake_prct_move(hh) = this_h5.prct_move*100;
    
    R_thresh = [R_thresh, this_h5.REM_pre_rThresh , this_h5.REM_post_rThresh];
    
    nPre_A(hh) = this_h5.pre_rem_nA;
    nPre_Act(hh) = sum(this_h5.pre_rem_nWake_A > 0);
    nPre_A_prct(hh) = (nPre_Act(hh)/nPre_A(hh))*100;
    nPre_sig(hh) = this_h5.pre_rem_A_sig_react;
    nPost_sig(hh) = this_h5.post_rem_A_sig_react;
    
    ReAct_str(hh) = this_h5.ReAct_Str_all_post;
    
    nPre_SWS_A(hh) = double(this_h5.pre_sws_A_sig_react);
    nPost_SWS_A(hh) = double(this_h5.post_sws_A_sig_react);
    
    
end

% table for stats

tbl = table(Sub', Cond', nWake_A', nPre_A', nPre_Act', nPre_A_prct', nPre_sig', nPost_sig', nPre_SWS_A', nPost_SWS_A', ReAct_str',nWake_R',wake_move', ...
    'VariableNames',{'Sub', 'Cond', 'nWake_A', 'nPre_A', 'nPre_Act', 'nPre_A_prct', 'nPre_sig', 'nPost_sig', 'nPre_SWS_A', 'nPost_SWS_A', 'ReAct_str', 'nWake_R', 'wake_move'});

tbl.Sub = nominal(tbl.Sub);
tbl.Cond = nominal(tbl.Cond);



%% plot number of significant assemblies
c_ord = MS_linspecer(14);
f_ord = [[253, 22, 26];[255, 179, 13]; [132, 147, 35];[67, 127, 151]; [0, 42, 95]]/255;

y_lim = [0 max([nWake_A, nPre_A])];



%%
% alt
f_ord = [[253, 22, 26];[255, 179, 13]; [132, 147, 35];[67, 127, 151]; [0, 42, 95]]/255;
f_post = [f_ord ,[.5, .5, .5, .5, .5]'];
off_set = 0:2:8;

figure(888)
clf
set(gcf, 'Units', 'centimeters', 'Position', [0 0  30 20])


% plot the pre and post REM assembly centroid reactivations.
subplot(2,4,5)
off_set = 0:2.5:10;
for ii = 1:5
    if ii == 1
        this_idx = lt1_idx;
    elseif ii == 2
        this_idx = lt5_idx;
    elseif ii == 3
        this_idx = H1_idx;
    elseif ii == 4
        this_idx = H5_idx;
    elseif ii == 5
        this_idx = HS_idx;
    end
    % add 'walls'
    if ii == 1 || ii == 2
        rectangle('position', [off_set(ii)-1, p_centr(1), 2, p_centr(end)- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii == 3 || ii ==4
        rectangle('position', [off_set(ii)-1, p_centr(ceil(length(p_centr)/2)), 2, p_centr(end)- p_centr(ceil(length(p_centr)/2))], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii ==5
        rectangle('position', [off_set(ii)-1, p_centr(1), 2, p_centr(ceil(length(p_centr)/2))- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    end
    
    hold on
    norm_val= max([nanmean(A_hist_pre(this_idx,:)) nanmean(A_hist_post(this_idx,:))]);
    plot((-nanmean(A_hist_pre(this_idx,:))./norm_val)+off_set(ii), p_centr,   'linewidth', 2, 'Color', f_ord(ii,:))
    plot((nanmean(A_hist_post(this_idx,:))./norm_val)+off_set(ii), p_centr,  'linewidth', 2, 'Color', [f_ord(ii,:) .5])
    
    if ii == 1
        %         text(off_set(ii)-.1, 97, '\leftarrow', 'HorizontalAlignment','right', 'Interpreter', 'TeX', 'fontsize', 8)
        text(off_set(ii)-.1, p_centr(end), 'pre', 'HorizontalAlignment','right', 'Interpreter', 'TeX', 'fontsize', 8, 'VerticalAlignment', 'baseline')
        %         text(off_set(ii)+.1, 97, '\rightarrow', 'HorizontalAlignment','left', 'Interpreter', 'TeX', 'fontsize', 8)
        text(off_set(ii)+.1, p_centr(end), 'post', 'HorizontalAlignment','left', 'Interpreter', 'TeX', 'fontsize', 8, 'VerticalAlignment', 'baseline')
        
    end
end
set(gca, 'xtick', off_set, 'ytick', [0 100])
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
ylabel('position (cm)')
axis('square')
title('REM assembly tuning')




subplot(2,4,1)
off_set = 0:2:8;
for ii = 1:5
    if ii == 1
        this_idx = lt1_idx;
    elseif ii == 2
        this_idx = lt5_idx;
    elseif ii == 3
        this_idx = H1_idx;
    elseif ii == 4
        this_idx = H5_idx;
    elseif ii == 5
        this_idx = HS_idx;
    end
    
    if ii == 1 || ii == 2
        rectangle('position', [off_set(ii), p_centr(1), mode(diff(off_set))*.75, p_centr(end)- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii == 3 || ii ==4
        rectangle('position', [off_set(ii), p_centr(ceil(length(p_centr)/2)), mode(diff(off_set))*.75, p_centr(end)- p_centr(ceil(length(p_centr)/2))], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii ==5
        rectangle('position', [off_set(ii), p_centr(1), mode(diff(off_set))*.75, p_centr(ceil(length(p_centr)/2))- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    end
    
    hold on
    norm_val= max([nanmean(A_hist_wake(this_idx,:))]);
    plot((nanmean(A_hist_wake(this_idx,:))./norm_val)+off_set(ii), p_centr,   'linewidth', 2, 'Color', f_ord(ii,:))
    
    
end
set(gca, 'xtick', off_set, 'ytick', [0 100])
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
ylabel('position (cm)')
axis('square')
title('Wake assembly tuning')


% print(gcf,  [fig_dir filesep   'Spatial_summery_' strrep(num2str(B_out{1}{1}.info.bin), '.', 'p') 's_bin.pdf'], '-dpdf','-r0')



%
% figure(99)
% clf
set(gcf, 'Units', 'centimeters', 'Position', [0 0  30 20])
subplot(2,4,2)
% MS_bar_w_err(nWake_A(Cond==1), nWake_A(Cond ==2), c_ord(3,:))
hb = bar([nanmean(nWake_A(Cond==1)), nanmean(nWake_A(Cond ==2)), nanmean(nWake_A(Cond ==3)), nanmean(nWake_A(Cond ==4)), nanmean(nWake_A(Cond ==5))], 'FaceColor', 'flat');
for ii = 1:5
    hb.CData(ii,:) = f_ord(ii,:);
end
hold on
eb = errorbar([nanmean(nWake_A(Cond==1)), nanmean(nWake_A(Cond ==2)), nanmean(nWake_A(Cond ==3)), nanmean(nWake_A(Cond ==4)), nanmean(nWake_A(Cond ==5))], [MS_SEM(nWake_A(Cond==1)), MS_SEM(nWake_A(Cond ==2)), MS_SEM(nWake_A(Cond ==3)), MS_SEM(nWake_A(Cond ==4)), MS_SEM(nWake_A(Cond ==5))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
ylim(y_lim);
title('Wake Assemblies in Wake')
axis('square')
ylabel('Number of assemblies')


subplot(2,4,6)
hb = bar([nanmean(nPre_A(Cond==1)), nanmean(nPre_A(Cond ==2)), nanmean(nPre_A(Cond ==3)), nanmean(nPre_A(Cond ==4)), nanmean(nPre_A(Cond ==5))] , 'FaceColor', 'flat');
for ii = 1:5
    hb.CData(ii,:) = f_ord(ii,:);
end
hold on

eb = errorbar([nanmean(nPre_A(Cond==1)), nanmean(nPre_A(Cond ==2)), nanmean(nPre_A(Cond ==3)), nanmean(nPre_A(Cond ==4)), nanmean(nPre_A(Cond ==5))], [MS_SEM(nPre_A(Cond==1)), MS_SEM(nPre_A(Cond ==2)), MS_SEM(nPre_A(Cond ==3)), MS_SEM(nPre_A(Cond ==4)), MS_SEM(nPre_A(Cond ==5))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
ylim(y_lim);
title('Pre Assemblies in Wake')
axis('square')

% subplot(1,4,2)
% % MS_bar_w_err(nPre_sig(Cond==1), nPre_sig(Cond ==2), c_ord(6,:))
% hb = bar([nanmean(nPre_sig(Cond==1)), nanmean(nPre_sig(Cond ==2)), nanmean(nPre_sig(Cond ==3)), nanmean(nPre_sig(Cond ==4)), nanmean(nPre_sig(Cond ==5))], 'FaceColor', 'flat');
% hold on
% for ii = 1:5
% hb.CData(ii,:) = f_ord(ii,:);
% end
% eb = errorbar([nanmean(nPre_sig(Cond==1)), nanmean(nPre_sig(Cond ==2)), nanmean(nPre_sig(Cond ==3)), nanmean(nPre_sig(Cond ==4)), nanmean(nPre_sig(Cond ==5))], [MS_SEM(nPre_sig(Cond==1)), MS_SEM(nPre_sig(Cond ==2)), MS_SEM(nPre_sig(Cond ==3)), MS_SEM(nPre_sig(Cond ==4)), MS_SEM(nPre_sig(Cond ==5))]);
% eb.LineStyle = 'none';
% eb.Color = 'k';
% set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
% ylim(y_lim);
% title('Wake Assemblies in Pre REM')
% axis('square')
%
%
% subplot(1,4,3)
%
% % MS_bar_w_err(nPost_sig(Cond==1), nPost_sig(Cond ==2), c_ord(2,:))
% hb = bar([nanmean(nPost_sig(Cond==1)), nanmean(nPost_sig(Cond ==2)), nanmean(nPost_sig(Cond ==3)), nanmean(nPost_sig(Cond ==4)), nanmean(nPost_sig(Cond ==5))], 'FaceColor', 'flat');
% for ii = 1:5
% hb.CData(ii,:) = f_ord(ii,:);
% end
% hold on
% eb = errorbar([nanmean(nPost_sig(Cond==1)), nanmean(nPost_sig(Cond ==2)), nanmean(nPost_sig(Cond ==3)), nanmean(nPost_sig(Cond ==4)), nanmean(nPost_sig(Cond ==5))], [MS_SEM(nPost_sig(Cond==1)), MS_SEM(nPost_sig(Cond ==2)), MS_SEM(nPost_sig(Cond ==3)), MS_SEM(nPost_sig(Cond ==4)), MS_SEM(nPost_sig(Cond ==5))]);
% eb.LineStyle = 'none';
% eb.Color = 'k';
% set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
% ylim(y_lim);
% axis('square')
% title('Wake Assemblies in Post REM')



% try to get the pre and post REM assemblies in one plot
pre_mean = [nanmean(nPre_sig(Cond==1)), nanmean(nPre_sig(Cond ==2)), nanmean(nPre_sig(Cond ==3)), nanmean(nPre_sig(Cond ==4)), nanmean(nPre_sig(Cond ==5))];
post_mean = [nanmean(nPost_sig(Cond==1)), nanmean(nPost_sig(Cond ==2)), nanmean(nPost_sig(Cond ==3)), nanmean(nPost_sig(Cond ==4)), nanmean(nPost_sig(Cond ==5))];
pre_err = [MS_SEM(nPre_sig(Cond==1)), MS_SEM(nPre_sig(Cond ==2)), MS_SEM(nPre_sig(Cond ==3)), MS_SEM(nPre_sig(Cond ==4)), MS_SEM(nPre_sig(Cond ==5))];
post_err = [MS_SEM(nPost_sig(Cond==1)), MS_SEM(nPost_sig(Cond ==2)), MS_SEM(nPost_sig(Cond ==3)), MS_SEM(nPost_sig(Cond ==4)), MS_SEM(nPost_sig(Cond ==5))];

y = [pre_mean; post_mean]';
err = [pre_err; post_err]';


subplot(2,4,7)
cla;
hb = bar(y, 1); % get the bar handles
hold on;
% xpos = ;
pause(.5)

errorbar(hb(1).XData + hb(1).XOffset, y(:,1), err(:,1), 'LineStyle', 'none', ...
    'Color', 'k');
pause(.5)
errorbar(hb(2).XData + hb(2).XOffset, y(:,2), err(:,2), 'LineStyle', 'none', ...
    'Color', 'k');
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
ylim(y_lim);
axis('square')
title('Wake assemblies in Pre and Post REM')
legend({'pre', 'post'}, 'box', 'off')

% REACT Str mean(mean Post - mean Pre)
% MS_bar_w_err(nWake_A(Cond==1), nWake_A(Cond ==2), c_ord(3,:))
subplot(2,4,8)
cla
hb = bar([nanmean(ReAct_str(Cond==1)), nanmean(ReAct_str(Cond ==2)), nanmean(ReAct_str(Cond ==3)), nanmean(ReAct_str(Cond ==4)), nanmean(ReAct_str(Cond ==5))], 'FaceColor', 'flat');
for ii = 1:5
    hb.CData(ii,:) = f_ord(ii,:);
end
hold on
eb = errorbar([nanmean(ReAct_str(Cond==1)), nanmean(ReAct_str(Cond ==2)), nanmean(ReAct_str(Cond ==3)), nanmean(ReAct_str(Cond ==4)), nanmean(ReAct_str(Cond ==5))],...
    [MS_SEM(ReAct_str(Cond==1)), MS_SEM(ReAct_str(Cond ==2)), MS_SEM(ReAct_str(Cond ==3)), MS_SEM(ReAct_str(Cond ==4)), MS_SEM(ReAct_str(Cond ==5))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
% ylim(y_lim);
axis('square')
title('Post REM Reactivation Strength')


% SetFigure([], gcf)

%  get the REM rate and space


% subplot(1,4,2)
% cla
% hb = bar([nanmean(wake_move(Cond==1)), nanmean(wake_move(Cond ==2)), nanmean(wake_move(Cond ==3)), nanmean(wake_move(Cond ==4)), nanmean(wake_move(Cond ==5))], 'FaceColor', 'flat');
% for ii = 1:5
% hb.CData(ii,:) = f_ord(ii,:);
% end
% hold on
% eb = errorbar([nanmean(wake_move(Cond==1)), nanmean(wake_move(Cond ==2)), nanmean(wake_move(Cond ==3)), nanmean(wake_move(Cond ==4)), nanmean(wake_move(Cond ==5))],...
%     [MS_SEM(wake_move(Cond==1)), MS_SEM(wake_move(Cond ==2)), MS_SEM(wake_move(Cond ==3)), MS_SEM(wake_move(Cond ==4)), MS_SEM(wake_move(Cond ==5))]);
% eb.LineStyle = 'none';
% eb.Color = 'k';
% set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
% % ylim(y_lim);
% axis('square')
% title('Total movement (m)')

% subplot(2,4,3)
% cla
% hb = bar([nanmean(wake_move(Cond==1)), nanmean(wake_move(Cond ==2)), nanmean(wake_move(Cond ==3)), nanmean(wake_move(Cond ==4)), nanmean(wake_move(Cond ==5))], 'FaceColor', 'flat');
% for ii = 1:5
% hb.CData(ii,:) = f_ord(ii,:);
% end
% hold on
% eb = errorbar([nanmean(wake_move(Cond==1)), nanmean(wake_move(Cond ==2)), nanmean(wake_move(Cond ==3)), nanmean(wake_move(Cond ==4)), nanmean(wake_move(Cond ==5))],...
%     [MS_SEM(wake_move(Cond==1)), MS_SEM(wake_move(Cond ==2)), MS_SEM(wake_move(Cond ==3)), MS_SEM(wake_move(Cond ==4)), MS_SEM(wake_move(Cond ==5))]);
% eb.LineStyle = 'none';
% eb.Color = 'k';
% set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
% % ylim(y_lim);
% axis('square')
% title('Total movement (m)')


subplot(2,4,4)
cla
hb = bar([nanmean(nWake_R(Cond==1)), nanmean(nWake_R(Cond ==2)), nanmean(nWake_R(Cond ==3)), nanmean(nWake_R(Cond ==4)), nanmean(nWake_R(Cond ==5))], 'FaceColor', 'flat');
for ii = 1:5
    hb.CData(ii,:) = f_ord(ii,:);
end
hold on
eb = errorbar([nanmean(nWake_R(Cond==1)), nanmean(nWake_R(Cond ==2)), nanmean(nWake_R(Cond ==3)), nanmean(nWake_R(Cond ==4)), nanmean(nWake_R(Cond ==5))],...
    [MS_SEM(nWake_R(Cond==1)), MS_SEM(nWake_R(Cond ==2)), MS_SEM(nWake_R(Cond ==3)), MS_SEM(nWake_R(Cond ==4)), MS_SEM(nWake_R(Cond ==5))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
% ylim(y_lim);
axis('square')
title('wake activation rate (Hz)')


subplot(2,4,3)
cla
hb = bar([nanmean(wake_prct_move(Cond==1)), nanmean(wake_prct_move(Cond ==2)), nanmean(wake_prct_move(Cond ==3)), nanmean(wake_prct_move(Cond ==4)), nanmean(wake_prct_move(Cond ==5))], 'FaceColor', 'flat');
for ii = 1:5
    hb.CData(ii,:) = f_ord(ii,:);
end
hold on
eb = errorbar([nanmean(wake_prct_move(Cond==1)), nanmean(wake_prct_move(Cond ==2)), nanmean(wake_prct_move(Cond ==3)), nanmean(wake_prct_move(Cond ==4)), nanmean(wake_prct_move(Cond ==5))],...
    [MS_SEM(wake_prct_move(Cond==1)), MS_SEM(wake_prct_move(Cond ==2)), MS_SEM(wake_prct_move(Cond ==3)), MS_SEM(wake_prct_move(Cond ==4)), MS_SEM(wake_prct_move(Cond ==5))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
% ylim(y_lim);
axis('square')
title('Percent of time moving')


set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(gcf,  [fig_dir filesep   'Stats_summery_' strrep(num2str(A_out{1}{1}.info.bin), '.', 'p') 's_bin.pdf'], '-dpdf','-r0')

%% stats
bin_idx = nearest_idx([40 60], A_out{1}{1}.map{1}.bins);

for jj = size(A_hist_pre,1):-1:1
    
    A_hist_pre_n = A_hist_pre(jj,:)./max(A_hist_pre(jj,:));
    A_hist_post_n = A_hist_post(jj,:)./max(A_hist_post(jj,:));
    
    if Cond(jj) ~= 5
        pre_open(jj) = nanmean(A_hist_pre_n(1:8));
        pre_closed(jj) = nanmean(A_hist_pre_n(12:end));
        
        wake_open(jj) = nanmean(A_hist_wake(jj,1:8));
        wake_closed(jj) = nanmean(A_hist_wake(jj,12:end));
        
        post_open(jj) = nanmean(A_hist_post_n(1:8));
        post_closed(jj) = nanmean(A_hist_post_n(12:end));
        
        
        pre_map_open(jj) = nanmean(Pre_ReAct_mean_map(jj,1:bin_idx(1)));
        pre_map_closed(jj) = nanmean(Pre_ReAct_mean_map(jj,bin_idx:end));
        
        wake_map_open(jj) = nanmean(Wake_ReAct_mean_map(jj,1:bin_idx(1)));
        wake_map_closed(jj) = nanmean(Wake_ReAct_mean_map(jj,bin_idx(2):end));
        
        post_map_open(jj) = nanmean(Post_ReAct_mean_map(jj,1:bin_idx(1)));
        post_map_closed(jj) = nanmean(Post_ReAct_mean_map(jj,bin_idx(2):end));
        
    elseif Cond(jj) ==5
        pre_closed(jj) = nanmean(A_hist_pre_n(1:8));
        pre_open(jj) = nanmean(A_hist_pre_n(12:end));
        
        wake_closed(jj) = nanmean(A_hist_wake(jj,1:8));
        wake_open(jj) = nanmean(A_hist_wake(jj,12:end));
        
        post_closed(jj) = nanmean(A_hist_post_n(1:8));
        post_open(jj) = nanmean(A_hist_post_n(12:end));
        
        pre_map_open(jj) = nanmean(Pre_ReAct_mean_map(jj,1:bin_idx(1)));
        pre_map_closed(jj) = nanmean(Pre_ReAct_mean_map(jj,bin_idx:end));
        
        wake_map_open(jj) = nanmean(Wake_ReAct_mean_map(jj,1:bin_idx(1)));
        wake_map_closed(jj) = nanmean(Wake_ReAct_mean_map(jj,bin_idx(2):end));
        
        post_map_open(jj) = nanmean(Post_ReAct_mean_map(jj,1:bin_idx(1)));
        post_map_closed(jj) = nanmean(Post_ReAct_mean_map(jj,bin_idx(2):end));
        
    end
    
    d_temp = A_hist_post_n./A_hist_pre_n;
    
    diff_open(jj) = nanmean(d_temp(1,1:8));
    diff_closed(jj) = nanmean(d_temp(12:end));
end
ReAct_tbl = table(Sub', Cond', pre_open', pre_closed', post_open', post_closed',pre_map_open', pre_map_closed',wake_map_open', wake_map_closed',post_map_open', post_map_closed', ...
    'VariableNames',{'Sub', 'Cond', 'pre_open', 'pre_closed', 'post_open', 'post_closed','pre_map_open', 'pre_map_closed','wake_map_open', 'wake_map_closed','post_map_open', 'post_map_closed'});

ReAct_tbl.Sub = nominal(ReAct_tbl.Sub);
ReAct_tbl.Cond = nominal(ReAct_tbl.Cond);


% long form.
all_hist_O_C = []; all_map_O_C = []; Pre_Post = {};  O_C = {}; S_l = []; S_type = [];

for jj = 1:length(ReAct_tbl.Sub)
    
    all_hist_O_C(end+1) = ReAct_tbl.pre_open(jj);
    all_hist_O_C(end+1) = ReAct_tbl.pre_closed(jj);
    all_hist_O_C(end+1) = ReAct_tbl.post_open(jj);
    all_hist_O_C(end+1) = ReAct_tbl.post_closed(jj);
    
    all_map_O_C(end+1) = ReAct_tbl.pre_map_open(jj);
    all_map_O_C(end+1) = ReAct_tbl.pre_map_closed(jj);
    all_map_O_C(end+1) = ReAct_tbl.post_map_open(jj);
    all_map_O_C(end+1) = ReAct_tbl.post_map_closed(jj);
    
    Pre_Post{end+1} = 'pre';
    Pre_Post{end+1} = 'pre';
    Pre_Post{end+1} = 'post';
    Pre_Post{end+1} = 'post';
    
    O_C{end+1} = 'open';
    O_C{end+1} = 'closed';
    O_C{end+1} = 'open';
    O_C{end+1} = 'closed';
    
    S_l(end+1) = ReAct_tbl.Sub(jj);
    S_l(end+1) = ReAct_tbl.Sub(jj);
    S_l(end+1) = ReAct_tbl.Sub(jj);
    S_l(end+1) = ReAct_tbl.Sub(jj);
    
    S_type(end+1) = ReAct_tbl.Cond(jj);
    S_type(end+1) = ReAct_tbl.Cond(jj);
    S_type(end+1) = ReAct_tbl.Cond(jj);
    S_type(end+1) = ReAct_tbl.Cond(jj);
    
end
H1_idx = S_type == 3;
L5_idx = S_type == 2;


ReAct_long_tbl = table(S_l(L5_idx |H1_idx)', S_type(L5_idx |H1_idx)',Pre_Post(L5_idx |H1_idx)', O_C(L5_idx |H1_idx)',all_map_O_C(L5_idx |H1_idx)', all_hist_O_C(L5_idx |H1_idx)', 'VariableNames', {'Subject', 'Session','Pre_post', 'Open_Closed', 'Mean_map', 'Mean_hist'});
ReAct_long_tbl.Subject = nominal(ReAct_long_tbl.Subject);
ReAct_long_tbl.Session = categorical(ReAct_long_tbl.Session);
ReAct_long_tbl.Open_Closed = categorical(ReAct_long_tbl.Open_Closed);
ReAct_long_tbl.Pre_post = categorical(ReAct_long_tbl.Pre_post);

LM_map = fitlme(ReAct_long_tbl, 'Mean_map ~ Session * Open_Closed * Pre_post + (1|Subject)')

figure(8989); clf
boxchart(ReAct_long_tbl.Session, ReAct_long_tbl.Mean_map, 'GroupbyColor', ReAct_long_tbl.Open_Closed)
legend
set(gca, 'xticklabel', {'Familiar', 'Anxiety'})

writetable(ReAct_long_tbl, [inter_dir filesep 'Map_ReAct_tbl.csv'])

% same but with the difference between open and closed
% long form.
all_hist_O_C = []; all_map_O_C = []; Pre_Post = {};  S_l = []; S_type = [];

for jj = 1:length(ReAct_tbl.Sub)
    
    all_hist_O_C(end+1) = ReAct_tbl.pre_open(jj) - ReAct_tbl.pre_closed(jj);
    all_hist_O_C(end+1) = ReAct_tbl.post_open(jj) - ReAct_tbl.post_closed(jj);
    
    all_map_O_C(end+1) = ReAct_tbl.pre_map_open(jj) - ReAct_tbl.pre_map_closed(jj);
    all_map_O_C(end+1) = ReAct_tbl.post_map_open(jj) - ReAct_tbl.post_map_closed(jj);
    
    Pre_Post{end+1} = 'pre';
    Pre_Post{end+1} = 'post';
    
    
    S_l(end+1) = ReAct_tbl.Sub(jj);
    S_l(end+1) = ReAct_tbl.Sub(jj);
    
    S_type(end+1) = ReAct_tbl.Cond(jj);
    S_type(end+1) = ReAct_tbl.Cond(jj);
    
end
H_idx = S_type == 3;
F_idx = S_type == 2;

% S_type(S_type == 2) = '


Diff_tbl = table(S_l(F_idx |H_idx)', S_type(F_idx |H_idx)',Pre_Post(F_idx |H_idx)', all_map_O_C(F_idx |H_idx)', all_hist_O_C(F_idx |H_idx)', 'VariableNames', {'Subject', 'Session','Pre_post', 'Mean_map', 'Mean_hist'});
Diff_tbl.Subject = nominal(Diff_tbl.Subject);
Diff_tbl.Session = categorical(Diff_tbl.Session);
Diff_tbl.Pre_post = categorical(Diff_tbl.Pre_post);


diff_LM_map = fitlme(Diff_tbl, 'Mean_map ~ Session * Pre_post + (1|Subject)')
anova(diff_LM_map)
figure(89892); clf
boxchart(Diff_tbl.Session, Diff_tbl.Mean_map, 'GroupbyColor', Diff_tbl.Pre_post)
legend
set(gca, 'xticklabel', {'Familiar', 'Anxiety'})
ylabel({'reactivation bias'; '<---- Closed    |    Open --->'})
writetable(Diff_tbl, [inter_dir filesep 'diff_tbl.csv'])



% rma_stats = ranova(
% rm = fitrm(ReAct_tbl, 'c2_t1-c1_t0 ~ 1', 'WithinDesign', w);
% ranova(rm, 'withinmodel', 'cond*time')
%% normalized reactivation maps


figure(8989)
clf
set(gcf, 'Units', 'centimeters', 'Position', [0 0  30 20])


% plot the pre and post REM assembly centroid reactivations.
off_set = 0:2.5:10;
for ii = 1:5
    if ii == 1
        this_idx = lt1_idx;
    elseif ii == 2
        this_idx = lt5_idx;
    elseif ii == 3
        this_idx = H1_idx;
    elseif ii == 4
        this_idx = H5_idx;
    elseif ii == 5
        this_idx = HS_idx;
    end
    % add 'walls'
    if ii == 1 || ii == 2
        rectangle('position', [off_set(ii)-1, p_centr(1), 2, p_centr(end)- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii == 3 || ii ==4
        rectangle('position', [off_set(ii)-1, p_centr(ceil(length(p_centr)/2)), 2, p_centr(end)- p_centr(ceil(length(p_centr)/2))], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii ==5
        rectangle('position', [off_set(ii)-1, p_centr(1), 2, p_centr(ceil(length(p_centr)/2))- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    end
    
    hold on
    %     norm_val= max([nanmean(A_hist_pre(this_idx,:)) nanmean(A_hist_post(this_idx,:))]);
    %     plot((-nanmean(A_hist_pre(this_idx,:))./norm_val)+off_set(ii), p_centr,   'linewidth', 2, 'Color', f_ord(ii,:))
    %     plot((nanmean(A_hist_post(this_idx,:))./norm_val)+off_set(ii), p_centr,  'linewidth', 2, 'Color', [f_ord(ii,:) .5])
    
    this_pre = A_hist_pre(this_idx,:);
    this_pre(this_pre == 0) = inf;
    this_post = A_hist_post(this_idx,:);
    this_post(this_post == 0) = inf;
    
    post_pre = mean(this_post./this_pre, 'omitnan');
    
    %     post_pre = nanmean(A_hist_post(this_idx,:)./A_hist_pre(this_idx,:));
    
    post_pre(isinf(post_pre)) = NaN;
    
    post_pre = post_pre./max(post_pre);
    
    post_pre(isnan(post_pre)) =0;
    
    plot(post_pre+off_set(ii), p_centr,'linewidth', 2, 'Color', f_ord(ii,:))
    xline(off_set(ii)+1)
end
set(gca, 'xtick', off_set, 'ytick', [0 100])
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
ylabel('position (cm)')
axis('square')
title('REM Post/Pre assembly tuning')
%%
% try to get the pre and post REM assemblies in one plot

figure(881)
pre_mean = [nanmean(nPre_sig(Cond==1)), nanmean(nPre_sig(Cond ==2)), nanmean(nPre_sig(Cond ==3)), nanmean(nPre_sig(Cond ==4)), nanmean(nPre_sig(Cond ==5))];
post_mean = [nanmean(nPost_sig(Cond==1)), nanmean(nPost_sig(Cond ==2)), nanmean(nPost_sig(Cond ==3)), nanmean(nPost_sig(Cond ==4)), nanmean(nPost_sig(Cond ==5))];
pre_err = [MS_SEM(nPre_sig(Cond==1)), MS_SEM(nPre_sig(Cond ==2)), MS_SEM(nPre_sig(Cond ==3)), MS_SEM(nPre_sig(Cond ==4)), MS_SEM(nPre_sig(Cond ==5))];
post_err = [MS_SEM(nPost_sig(Cond==1)), MS_SEM(nPost_sig(Cond ==2)), MS_SEM(nPost_sig(Cond ==3)), MS_SEM(nPost_sig(Cond ==4)), MS_SEM(nPost_sig(Cond ==5))];

y = [pre_mean; post_mean]';
err = [pre_err; post_err]';
subplot(1,4,3)
cla;

hb = bar(y, 1); % get the bar handles
hold on;

% xpos = ;
pause(.5)

errorbar(hb(1).XData + hb(1).XOffset, y(:,1), err(:,1), 'LineStyle', 'none', ...
    'Color', 'k');
pause(.5)
errorbar(hb(2).XData + hb(2).XOffset, y(:,2), err(:,2), 'LineStyle', 'none', ...
    'Color', 'k');
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
ylim(y_lim);
axis('square')
title('Wake assemblies in Pre and Post REM')

%% quick stats


disp('Nov FAM AnxN AnxF AnxS')
lme_mod = fitlme(tbl,'nWake_A~Cond+(1|Sub)', 'DummyVarCoding', 'reference');
A = lme_mod.anova; 
fprintf('%s |  F: %2.4f test p = %2.4f\n',lme_mod.Formula,   A.FStat(2), coefTest(lme_mod))
% 
% lme_mod = fitlme(tbl,'nPre_sig~Cond+(1|Sub)');
% A = lme_mod.anova; 
% fprintf('%s |  F: %2.4f test p = %2.4f\n',lme_mod.Formula,   A.FStat(2), coefTest(lme_mod))
% 
% lme_mod = fitlme(tbl,'nPost_sig~Cond+(1|Sub)');
% A = lme_mod.anova; 
% fprintf('%s |  F: %2.4f test p = %2.4f\n',lme_mod.Formula,   A.FStat(2), coefTest(lme_mod))
% 
% lme_mod = fitlme(tbl,'nPre_A~Cond+(1|Sub)');
% A = lme_mod.anova; 
% fprintf('%s |  F: %2.4f test p = %2.4f\n',lme_mod.Formula,   A.FStat(2), coefTest(lme_mod))

lme_mod = fitlme(tbl,'ReAct_str ~Cond+(1|Sub)', 'DummyVarCoding', 'reference');
A = lme_mod.anova; 
fprintf('%s |  F: %2.4f test p = %2.4f\n',lme_mod.Formula,   A.FStat(2), coefTest(lme_mod))

lme_mod = fitlme(tbl,'nWake_R ~Cond+(1|Sub)', 'DummyVarCoding', 'reference');
A = lme_mod.anova; 
fprintf('%s |  F: %2.4f test p = %2.4f\n',lme_mod.Formula,   A.FStat(2), coefTest(lme_mod))

NFA_tbl = tbl; 

NFA_tbl((NFA_tbl.Cond == '4' | NFA_tbl.Cond == '5'), :) = []; 
NFA_tbl.Cond = nominal(NFA_tbl.Cond);

disp('Novel, Fam, Anxiety only')
lme_mod = fitlme(NFA_tbl,'nWake_A~Cond+(1|Sub)', 'DummyVarCoding', 'reference');
A = lme_mod.anova; 
fprintf('%s |  F: %2.4f test p = %2.4f\n',lme_mod.Formula,   A.FStat(2), coefTest(lme_mod))
% 
% lme_mod = fitlme(tbl,'nPre_sig~Cond+(1|Sub)');
% A = lme_mod.anova; 
% fprintf('%s |  F: %2.4f test p = %2.4f\n',lme_mod.Formula,   A.FStat(2), coefTest(lme_mod))
% 
% lme_mod = fitlme(tbl,'nPost_sig~Cond+(1|Sub)');
% A = lme_mod.anova; 
% fprintf('%s |  F: %2.4f test p = %2.4f\n',lme_mod.Formula,   A.FStat(2), coefTest(lme_mod))
% 
% lme_mod = fitlme(tbl,'nPre_A~Cond+(1|Sub)');
% A = lme_mod.anova; 
% fprintf('%s |  F: %2.4f test p = %2.4f\n',lme_mod.Formula,   A.FStat(2), coefTest(lme_mod))

lme_mod = fitlme(NFA_tbl,'ReAct_str ~Cond+(1|Sub)', 'DummyVarCoding', 'reference');
A = lme_mod.anova; 
fprintf('%s |  F: %2.4f test p = %2.4f\n',lme_mod.Formula,   A.FStat(2), coefTest(lme_mod))

lme_mod = fitlme(NFA_tbl,'nWake_R ~Cond+(1|Sub)', 'DummyVarCoding', 'reference');
A = lme_mod.anova; 
fprintf('%s |  F: %2.4f test p = %2.4f\n',lme_mod.Formula,   A.FStat(2), coefTest(lme_mod))
% stats = anova1(lme_wake_A);
%  multcompare(stats)

% anova2(tbl, {nPre_sig, nPost_sig}, 'Cond')

% same but with ANOVAs



%% GLM comparing Novel and limit to just familiar vs Anx

tbl_glm = tbl; 


lt_idx = (tbl.Cond == '2'); 
ht_idx = tbl.Cond == '3'; 

tbl_lt_ht = tbl; 
tbl_lt_ht(~(lt_idx | ht_idx),:) = []; 


lt_ht_nPre_A = fitlme(tbl,'nPre_A  ~ Cond + (Sub)');
fprintf('LME - nPre_A: coef test p = %2.4f\n', coefTest(lt_ht_nPre_A))


% new var for LT5 vs HAT1

%% check against anovan
y = [tbl_lt_ht.nPre_sig ; tbl_lt_ht.nPost_sig];
cond = ones(size(tbl_lt_ht.Cond)); 
cond(tbl_lt_ht.Cond == '2') = 1; 
cond(tbl_lt_ht.Cond == '3') = 2; 
cond = [cond; cond]; 

cond_c = cell(size(cond)); 
for ii = length(cond_c)
    if cond(ii) == 1
        cond_c{ii} = 'Fam'; 
    elseif cond(ii) == 2;
        cond_c{ii} = 'Anx';
    end
end

pre_post = [ones(size(tbl_lt_ht.Cond));  ones(size(tbl_lt_ht.Cond))+1];
% g1 = tbl_lt_ht.Cond; 

% % tbl_lt_ht = [];
% tbl_

% lme_pre_post_sig = fitlme(tbl_lt_ht, 'score'

[p, tbl, stats] = anovan(y, {cond,pre_post}, 2,2, {'cond', 'pre_post'});

c1 = multcompare(stats);
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])



figure(1919)
subplot(2,2,1)
boxchart(pre_post, y,'GroupByColor',cond)
subplot(2,2,3:4)
bar([median(y(cond == 1 & pre_post == 1)), median(y(cond == 1 & pre_post == 2)); median(y(cond == 2 & pre_post == 1)), median(y(cond == 2 & pre_post == 2))]) 
hold on
errorbar_groups([median(y(cond == 1 & pre_post == 1)), median(y(cond == 1 & pre_post == 2)); median(y(cond == 2 & pre_post == 1)), median(y(cond == 2 & pre_post == 2))])
% anova2(tbl.nWake_A, tbl_lt_ht.Cond)


%% plot number of significant assemblies

y_lim = [0 100];



subplot(2,4,5)
MS_bar_w_err(nPre_A(Cond==1)./nWake_A(Cond==1)*100, nPre_A(Cond ==2)./nWake_A(Cond==2)*100, c_ord(12,:))
set(gca, 'XTickLabel', {'Novel', 'Familiar'})
ylim(y_lim);
title('Pre Assemblies in Wake')

ylabel('% of wake assemblies')

subplot(2,4,7)
MS_bar_w_err(nWake_A(Cond==1)./nWake_A(Cond==1)*100, nWake_A(Cond ==2)./nWake_A(Cond==2)*100, c_ord(3,:))
set(gca, 'XTickLabel', {'Novel', 'Familiar'})
ylim(y_lim);

title('Wake Assemblies in Wake')


subplot(2,4,6)
MS_bar_w_err(nPre_sig(Cond==1)./nWake_A(Cond==1)*100, nPre_sig(Cond ==2)./nWake_A(Cond==2)*100, c_ord(6,:))
set(gca, 'XTickLabel', {'Novel', 'Familiar'})
ylim(y_lim);

title('Wake Assemblies in Pre REM')

subplot(2,4,8)
MS_bar_w_err(nPost_sig(Cond==1)./nWake_A(Cond==1)*100, nPost_sig(Cond ==2)./nWake_A(Cond==2)*100, c_ord(2,:))
set(gca, 'XTickLabel', {'Novel', 'Familiar'})
ylim(y_lim);

title('Wake Assemblies in Post REM')

cfg_set= [];
cfg_set.resize = 0;
cfg_set.ft_size = 10;
SetFigure(cfg_set, gcf)

%% SWS

c_idx = (Cond ==2 | Cond == 1);

y_lim = [0 max([nPre_sig(c_idx), nPost_sig(c_idx), nPre_SWS_A(c_idx), nPost_SWS_A(c_idx)])];


figure(99)
clf
set(gcf, 'Units', 'centimeters', 'Position', [0 0  30 20])
subplot(2,2,1)
% hb = bar([nanmean(nPre_A(Cond==1)),nanmean(nPre_sig(Cond==1)), nanmean(nPost_sig(Cond==1)), nanmean(nPre_SWS_A(Cond==1)), nanmean(nPost_SWS_A(Cond==1)); nanmean(nPre_A(Cond==2)),nanmean(nPre_sig(Cond==2)), nanmean(nPost_sig(Cond==2)), nanmean(nPre_SWS_A(Cond==2)), nanmean(nPost_SWS_A(Cond==1))]);
% hold on
% eb = errorbar([nanmean(nPre_A(Cond==1)),nanmean(nPre_sig(Cond==1)), nanmean(nPost_sig(Cond==1)), nanmean(nPre_SWS_A(Cond==1)), nanmean(nPost_SWS_A(Cond==1)); nanmean(nPre_A(Cond==2)),nanmean(nPre_sig(Cond==2)), nanmean(nPost_sig(Cond==2)), nanmean(nPre_SWS_A(Cond==2)), nanmean(nPost_SWS_A(Cond==1))],...
%     [MS_SEM(nPre_A(Cond==1)),MS_SEM(nPre_sig(Cond==1)), MS_SEM(nPost_sig(Cond==1)), MS_SEM(nPre_SWS_A(Cond==1)), MS_SEM(nPost_SWS_A(Cond==1)); MS_SEM(nPre_A(Cond==2)),MS_SEM(nPre_sig(Cond==2)), MS_SEM(nPost_sig(Cond==2)), MS_SEM(nPre_SWS_A(Cond==2)), MS_SEM(nPost_SWS_A(Cond==1))]);
% eb.LineStyle = 'none';
%
% set(gca, 'XTickLabel', {'Novel', 'Familiar'})
% legend({'Pre A in wake', 'Pre REM', 'Post REM', 'Pre SWS', 'Post SWS'}, 'Box', 'off', 'orientation', 'horizontal')
% hb(1).FaceColor = c_ord(12,:);
% hb(2).FaceColor = c_ord(6,:);
% hb(3).FaceColor = c_ord(3,:);
% hb(4).FaceColor = c_ord(9,:);
% hb(5).FaceColor = c_ord(10,:);



MS_bar_w_err(nPre_sig(Cond==1), nPre_sig(Cond ==2), c_ord(12,:))
set(gca, 'XTickLabel', {'Novel', 'Familiar'})
ylim(y_lim);
title('Wake assemblies in Pre REM')

ylabel('Number of assemblies')

subplot(2,2,3)

MS_bar_w_err(nPost_sig(Cond==1), nPost_sig(Cond ==2), c_ord(10,:))
set(gca, 'XTickLabel', {'Novel', 'Familiar'})
ylim(y_lim);

title('Wake Assemblies in Post REM')



subplot(2,2,2)
MS_bar_w_err(nPre_SWS_A(Cond==1), nPre_SWS_A(Cond ==2), c_ord(3,:))
set(gca, 'XTickLabel', {'Novel', 'Familiar'})
ylim(y_lim);

title('Wake assemblies in Pre SWS')

subplot(2,2,4)
MS_bar_w_err(nPost_SWS_A(Cond==1), nPost_SWS_A(Cond ==2), c_ord(4,:))
set(gca, 'XTickLabel', {'Novel', 'Familiar'})
ylim(y_lim);

title('Wake Assemblies in Post SWS')


%% plot a specific example
N_idx = []; F_idx = [];
for ii = length(A_out):-1:1
    
    if strcmpi(A_out{ii}{1}.info.subject, 'pv1060') &&  strcmpi(A_out{ii}{1}.info.session, 'LTD1')
        N_idx(ii) = true;
        F_idx(ii) = false;
    elseif strcmpi(A_out{ii}{1}.info.subject, 'pv1060') &&  strcmpi(A_out{ii}{1}.info.session, 'LTD5')
        N_idx(ii) = false;
        F_idx(ii) = true;
    else
        N_idx(ii) = false;
        F_idx(ii) = false;
    end
    
end

F_idx = find(F_idx);
N_idx = find(N_idx);


figure(9909)
clf
for iD = 1:2
    
    if iD == 1
        this_data = A_out{N_idx}{1};
    elseif iD ==2
        this_data = A_out{F_idx}{1};
    end
    
    p_cent = []; p_rate = [];
    
    for ii = length(this_data.P_pos):-1:1
        p_cent(ii) = this_data.map{ii}.cent_z;
        p_rate(ii) = length(this_data.P_loc{ii}.loc_time);
    end
    
    p_rate(p_cent > -1.95) = 0;
    [~, s_idx] = sort(p_rate, 'descend');
    p_rank = s_idx;
    
    if length(this_data.P_pos) >3
        p_idx = p_rank(1:4);
    else
        p_idx = p_rank(1:length(this_data.P_pos));
    end
    
    c_ord = [67, 127, 151; 132 147 35; 255 179 13; 253 22 26]/255;
    
    ax(iD) = subplot(2,2,iD);
    %         yline(this_data.REM_Post_stats.R_thresh, '--', 'color', [.7 .7 .7], 'linewidth', 0.3)
    yline(log10(this_data.REM_Post_stats.R_thresh), '--', 'color', [.7 .7 .7], 'linewidth', 0.3)
    
    hold on
    for ii = length(p_idx):-1:1
        this_proj = this_data.REM_Post_proj(p_idx(ii),:);
        plot(this_data.REM_Post_tvec, log10(this_proj), 'color', c_ord(ii,:), 'linewidth', .5)
    end
    
    ylim([0 inf])
    y_val = get(gca, 'YTick');
    %         y_val = [0 1 2 2.3979];
    % yTick', [0 1 2 2.3979]
    set(gca, 'YTickLabel', 10.^y_val, 'linewidth', 1);
    
    if iD == 1
        ylabel({'assembly strength'})
    end
    set(gca, 'xtick', []);
    
    if iD == 1
        xlim([0 60]);
    elseif iD == 2
        %                 set(gca, 'xtick', 0:60:this_data.REM_Post_tvec(end))
        %         set(gca, 'xticklabel', get(gca, 'xtick') - 90)
        xlim([90 150]);
    end
    
    ax2(iD) = subplot(2,2,iD+2);
    for ii = length(p_idx):-1:1
        this_proj_idx = find(this_data.REM_Post_proj(p_idx(ii),:) > this_data.REM_Post_stats.R_thresh);
        line([this_data.REM_Post_tvec(this_proj_idx); this_data.REM_Post_tvec(this_proj_idx)] , [ii-.5; ii+.5], 'color', c_ord(ii,:), 'linewidth', .5)
    end
    
    set(gca, 'ytick', [])
    
    if iD == 1
        xlim([0 60]);
        xlabel('time (s)')
        set(gca, 'xtick', [0 60])
        
    elseif iD == 2
        %         set(gca, 'xtick', 0:60:this_data.REM_Post_tvec(end))
        %         set(gca, 'xticklabel', get(gca, 'xtick') - 90)
        xlim([90 150]);
        %         x_lim = get(gca, 'xlim');
        set(gca, 'xtick', [])
        %                             set(gca, 'xtick', [0 60])
        
    end
    
end
linkaxes(ax, 'y')
set(gcf,'PaperUnits','inches', 'Units', 'inches');
set(gcf, 'position', [5 5 7 2])
set(gcf,'PaperSize', [7, 2]);
set(gca,'xlimmode','manual','ylimmode','manual')
print(gcf, '-dpdf', [fig_dir filesep 'Fig3_example.pdf'])

% save an h5 with the data for plotting.
fname = ['Fig3_assembly_.h5'];

if exist(fname, 'file')
    delete(fname)
end

hdf5write(fname, '/mouse', string(A_out{F_idx}{1}.info.subject));
hdf5write(fname, '/N_condition', string(A_out{N_idx}{1}.info.session),'WriteMode', 'append');
hdf5write(fname, '/F_condition', string(A_out{F_idx}{1}.info.session),'WriteMode', 'append');

hdf5write(fname, '/N_REM_post_rThresh', A_out{N_idx}{1}.REM_Post_stats.R_thresh,'WriteMode', 'append');
hdf5write(fname, '/F_REM_post_rThresh', A_out{F_idx}{1}.REM_Post_stats.R_thresh,'WriteMode', 'append');
hdf5write(fname, '/N_REM_post_proj', A_out{N_idx}{1}.REM_Post_proj(p_idx,:),'WriteMode', 'append');
hdf5write(fname, '/F_REM_post_proj', A_out{F_idx}{1}.REM_Post_proj(p_idx,:),'WriteMode', 'append');
hdf5write(fname, '/N_REM_post_tvec', A_out{N_idx}{1}.REM_Post_tvec,'WriteMode', 'append');
hdf5write(fname, '/F_REM_post_tvec', A_out{F_idx}{1}.REM_Post_tvec,'WriteMode', 'append');

hdf5write(fname, '/N_REM_post_xlim', int16([90 150]),'WriteMode', 'append');
hdf5write(fname, '/F_REM_post_xlim', int16([0 60]),'WriteMode', 'append');



%% Fig4 'Pre REM' example

P_idx = [];
for ii = length(A_out):-1:1
    
    if strcmpi(A_out{ii}{1}.info.subject, 'pv1060') &&  strcmpi(A_out{ii}{1}.info.session, 'LTD1')
        P_idx(ii) = true;
    else
        P_idx(ii) = false;
    end
    
end

P_idx = find(P_idx);





this_data = A_out{P_idx}{1};


MS_Asmbly_plot_raster_ReAct(this_data,[], 'pREM_data',[1 2 5 6])


xlim([1 90])
set(gcf,'PaperUnits','inches', 'Units', 'inches');
set(gcf, 'position', [5 5 7 2])
set(gcf,'PaperSize', [7, 2]);
set(gca,'xlimmode','manual','ylimmode','manual')
print(gcf, '-dpdf', [fig_dir filesep 'Fig4_example.pdf'])

% post ASBMLY in wake

MS_Asmbly_plot_REM_wake_raster_figure(this_data,'Post',[])

%% ReActivation stats across training/testing sets
nReAct = NaN(3,3,length(A_out));
nA_ReAct = nReAct;
nReAct_Rate = nReAct;
A_ReAct_Rate = nReAct;

d_list = {'LTD1','LTD5', 'HATD1', 'HATD5', 'HATDSwitch'};
sub_list = {'pv1043', 'pv1060', 'pv1069', 'pv1191', 'pv1192', 'pv1252'};
for iA = length(A_out):-1:1
    
    D_idx = find(contains(d_list,A_out{iA}{1}.info.session));
    %     fprintf('%s\n', d_list{D_idx});
    S_idx = find(contains(sub_list,A_out{iA}{1}.info.subject));
    %     fprintf('%s - %s\n', A_out{iA}{1}.info.subject, sub_list{S_idx});
    this_ReAct = MS_Asmbly_ReAct_Matrix(A_out{iA}{1});
    
    
    nReAct(:,:,iA) = this_ReAct.nReAct./max(this_ReAct.nReAct, [], 'all');
    nA_ReAct(:,:,iA) = this_ReAct.nA_ReAct./max(this_ReAct.nA_ReAct, [], 'all');
    nReAct_Rate(:,:,iA) = this_ReAct.nReAct_Rate./max(this_ReAct.nReAct_Rate, [], 'all');
    A_ReAct_Rate(:,:,iA) = this_ReAct.A_ReAct_Rate./max(this_ReAct.A_ReAct_Rate, [], 'all');
    labels = this_ReAct.info.labels;
    
    
    
end

figure(7778)
m = 4;
subplot(m,4,1)
imagesc(nanmean(nReAct(:,:,(novel_idx == 1 & anx_idx == 0)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
title('number of reactivations'); colorbar;caxis([0 1]);
ylabel({'LTD1'; 'Training'});

subplot(m,4,2)
imagesc(nanmean(nA_ReAct(:,:,(novel_idx == 1 & anx_idx == 0)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
title('number of assemblies reactivated'); colorbar;caxis([0 1]);

subplot(m,4,3)
imagesc(nanmean(nReAct_Rate(:,:,(novel_idx == 1 & anx_idx == 0)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
title('Rate(nA ./ length of recording)');colorbar; caxis([0 1]);

subplot(m,4,4)
imagesc(nanmean(A_ReAct_Rate(:,:,(novel_idx == 1 & anx_idx == 0)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
title('ReAct Rate per assembly');     colorbar; caxis([0 1]);


% LTD5
subplot(m,4,5)
imagesc(nanmean(nReAct(:,:,(novel_idx == 0 & anx_idx == 0)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
colorbar; caxis([0 1]);
ylabel({'LTD5'; 'Training'});


subplot(m,4,6)
imagesc(nanmean(nA_ReAct(:,:,(novel_idx == 0 & anx_idx == 0)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
colorbar; caxis([0 1]);

subplot(m,4,7)
imagesc(nanmean(nReAct_Rate(:,:,(novel_idx == 0 & anx_idx == 0)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
colorbar; caxis([0 1]);

subplot(m,4,8)
imagesc(nanmean(A_ReAct_Rate(:,:,(novel_idx == 0 & anx_idx == 0)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
colorbar; caxis([0 1]);


% HAT1
subplot(m,4,9)

imagesc(nanmean(nReAct(:,:,(novel_idx == 1 & anx_idx == 1)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
colorbar; caxis([0 1]);
ylabel({'HAT1'; 'Training'});

subplot(m,4,10)
imagesc(nanmean(nA_ReAct(:,:,(novel_idx == 1 & anx_idx == 1)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
colorbar; caxis([0 1]);

subplot(m,4,11)
imagesc(nanmean(nReAct_Rate(:,:,(novel_idx == 1 & anx_idx == 1)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
colorbar; caxis([0 1]);

subplot(m,4,12)
imagesc(nanmean(A_ReAct_Rate(:,:,(novel_idx == 1 & anx_idx == 1)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
colorbar; caxis([0 1]);

% HAT5
subplot(m,4,13)

imagesc(nanmean(nReAct(:,:,(novel_idx == 0 & anx_idx == 1)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
colorbar; caxis([0 1]);
ylabel({'HAT5'; 'Training'});

subplot(m,4,14)
imagesc(nanmean(nA_ReAct(:,:,(novel_idx == 0 & anx_idx == 1)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
colorbar; caxis([0 1]);

subplot(m,4,15)
imagesc(nanmean(nReAct_Rate(:,:,(novel_idx == 0 & anx_idx == 1)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
colorbar; caxis([0 1]);

subplot(m,4,16)
imagesc(nanmean(A_ReAct_Rate(:,:,(novel_idx == 0 & anx_idx == 1)), 3));
set(gca, 'xtick', 1:3, 'xticklabels', {'post', 'wake', 'pre'}, 'ytick', 1:3, 'yticklabels', {'post', 'wake', 'pre'})
colorbar; caxis([0 1]);

%% check the bias in assemblies


for iA = length(A_out):-1:1
    
    for ii = size(A_out{iA}{1}.P_proj,1):-1:1
        
        this_map = A_out{iA}{1}.Place_map{ii};
        
%         if this_map.cent_z < -1.95
            %             p = polyshape(this_map.bins, this_map.map_mean);
            %             cent = centroid(p);
            
            
            c_mean = mean(this_map.bins(this_map.peak)*this_map.bin_size);
            
            c_group = this_map.bins(this_map.peak)*this_map.bin_size;
            
            % get the acitvity of each cell in the assembly.
            c_idx = A_out{iA}{1}.P_pos{ii};
            
            % check the pre activity per cell per reactivation.
            % Pre
            [~, R_idx] = findpeaks(A_out{iA}{1}.REM_Pre_proj(ii,:),'MinPeakHeight', A_out{iA}{1}.REM_Pre_stats.R_thresh);
            
            if ~isempty(R_idx)
                A_out{iA}{1}.REM_Pre_act.A_act_p{ii} = sum(A_out{iA}{1}.REM_Pre_data(R_idx,c_idx)>0,1) ./ length(R_idx);
                A_out{iA}{1}.REM_Pre_act.R_idx{ii} = R_idx;
                AxC_map = [];
                for iC = 1:length(c_idx)
                    AxC_map = [AxC_map; repmat(this_map.map(iC,:), A_out{iA}{1}.REM_Pre_act.A_act_p{ii}(iC)*length(A_out{iA}{1}.REM_Pre_act.R_idx{ii}),1)];
                end
                A_out{iA}{1}.REM_Pre_act.AxCmap{ii} = AxC_map;
                A_out{iA}{1}.REM_Pre_act.AxCmap_mean{ii} = mean(AxC_map);
            else
                 A_out{iA}{1}.REM_Pre_act.A_act_p{ii} = []; 
                 A_out{iA}{1}.REM_Pre_act.R_idx{ii} = []; 
                 A_out{iA}{1}.REM_Pre_act.AxCmap{ii} = []; 
                 A_out{iA}{1}.REM_Pre_act.AxCmap_mean{ii} = []; 
            end
            % for plotting checks.
            %            [~, c_sort] = sort(A_out{iA}{1}.P_temp(c_idx,ii), 'descend')
            %            A_act_mat = A_out{iA}{1}.REM_Pre_data(:,c_idx(c_sort));
            
            % wake
            [~, R_idx] = findpeaks(A_out{iA}{1}.P_proj(ii,:),'MinPeakHeight', 10);
            if ~isempty(R_idx)
                A_out{iA}{1}.Wake_act.A_act_p{ii} =  sum(A_out{iA}{1}.wake_data(R_idx,c_idx)>0,1) ./ length(R_idx);
                A_out{iA}{1}.Wake_act.R_idx{ii} = R_idx;
                
                AxC_map = [];
                for iC = 1:length(c_idx)
                    AxC_map = [AxC_map; repmat(this_map.map(iC,:), A_out{iA}{1}.Wake_act.A_act_p{ii}(iC)*length(A_out{iA}{1}.Wake_act.R_idx{ii}),1)];
                end
                A_out{iA}{1}.Wake_act.AxCmap{ii} = AxC_map;
                A_out{iA}{1}.Wake_act.AxCmap_mean{ii} = mean(AxC_map);
            else
                A_out{iA}{1}.Wake_act.A_act_p{ii} = [];
                A_out{iA}{1}.Wake_act.R_idx{ii} = [];
                A_out{iA}{1}.Wake_act.AxCmap{ii} = [];
                A_out{iA}{1}.Wake_act.AxCmap_mean{ii} = [];
            end
            

            % REM_post
            [~, R_idx] = findpeaks(A_out{iA}{1}.REM_Post_proj(ii,:),'MinPeakHeight', A_out{iA}{1}.REM_Post_stats.R_thresh);
            if ~isempty(R_idx)
                A_out{iA}{1}.REM_Post_act.A_act_p{ii} = sum(A_out{iA}{1}.REM_Post_data(R_idx,c_idx)>0,1) ./ length(R_idx);
                A_out{iA}{1}.REM_Post_act.R_idx{ii} = R_idx;
                
                AxC_map = [];
                for iC = 1:length(c_idx)
                    AxC_map = [AxC_map; repmat(this_map.map(iC,:), A_out{iA}{1}.REM_Post_act.A_act_p{ii}(iC)*length(A_out{iA}{1}.REM_Post_act.R_idx{ii}),1)];
                end
                A_out{iA}{1}.REM_Post_act.AxCmap{ii} = AxC_map;
                A_out{iA}{1}.REM_Post_act.AxCmap_mean{ii} = mean(AxC_map);
            else
                A_out{iA}{1}.REM_Post_act.A_act_p{ii} = [];
                A_out{iA}{1}.REM_Post_act.R_idx{ii} = [];
                A_out{iA}{1}.REM_Post_act.AxCmap{ii} = [];
                A_out{iA}{1}.REM_Post_act.AxCmap_mean{ii} = [];
            end
            
%         end
    end
end


%% get the range of REM R thresholds
pre_R = []; 
post_R = []; 
pREM_R = []; 

for iA = length(A_out):-1:1
    
        
        pre_R(end+1) = A_out{iA}{1}.REM_Pre_stats.R_thresh;
        
        
        post_R(end+1) = A_out{iA}{1}.REM_Post_stats.R_thresh;

        if isfield( A_out{iA}{1}, 'pREM_stats')
        pREM_R(end+1) = A_out{iA}{1}.pREM_stats.R_thresh; 
        end
end

fprintf(['Pre R Thresh: %0.2f (' char(177) '%0.2f)\n'], mean(pre_R), std(pre_R))
fprintf(['Post R Thresh: %0.2f (' char(177) '%0.2f)\n'], mean(post_R), std(post_R))
fprintf(['Pre & post R Thresh: %0.2f (' char(177) '%0.2f)\n'], mean([pre_R post_R]), std([pre_R post_R]))
fprintf(['All R Thresh: %0.2f (' char(177) '%0.2f)\n'], mean(pREM_R), std([pre_R post_R]))

% now what? Quantify the spatial bias? Plot things to see what might be
% changing? 
%% %%%%%%%%%%%%  sample plots for PCA ICA methods. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1111)

subplot(2,4,1)
cla

hold on

imagesc(zscore(A_out{ii}{1}.wake_data', [], 2))
ylim([1 256]); xlim([0 120])
axis('square')
%%
figure(919)
w_ord = MS_linspecer(5);

for ii = 1:5
    subplot(1,5,ii)
    cla
    % stem(, A_out{ii}{1}.P_temp(:,ii))
    hold on
    stem(1:256, A_out{ii}{1}.A_temp(:,ii+20), 'color', [.8 .8 .8 .2])
    % stem(A_out{ii}{1}.P_pos{ii}, A_out{ii}{1}.P_temp(A_out{ii}{1}.P_pos{ii},ii), 'color', w_ord(ii,:), 'MarkerFaceColor', w_ord(ii,:))
    view(90,90)
    xlim([1 256])
    ylim([-.5 0.7])
    set(gca, 'xticklabel', [])
    axis off
end

exportgraphics(gcf,[fig_dir filesep  'ica_null.pdf'],'ContentType','vector')


%% SFN poster plots


%% Specific example for Anxiety

%% bar graph of reactivations

% f_ord = [[253, 22, 26];[255, 179, 13]; [132, 147, 35];; [0, 42, 95]]/255;
f_ord = [0.2627 0.4980 0.5922; .7 .7 .7; 0.9922 0.0863 0.1020;1.0000 0.7020 0.0510; 0.5176 0.5765 0.1373]; 
f_post = [f_ord ,[.5, .5, .5, .5, .5]'];
off_set = 0:2:8;

figure(8118)
clf
set(gcf, 'Units', 'centimeters', 'Position', [0 0  30 20])


% plot the pre and post REM assembly centroid reactivations.
subplot(2,4,5)
off_set = 0:2.5:10;
for ii = 1:5
    if ii == 1
        this_idx = lt1_idx;
    elseif ii == 2
        this_idx = lt5_idx;
    elseif ii == 3
        this_idx = H1_idx;
    elseif ii == 4
        this_idx = H5_idx;
    elseif ii == 5
        this_idx = HS_idx;
    end
    % add 'walls'
    if ii == 1 || ii == 2
        rectangle('position', [off_set(ii)-1, p_centr(1), 2, p_centr(end)- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii == 5
        rectangle('position', [off_set(ii)-1, p_centr(ceil(length(p_centr)/2)), 2, p_centr(end)- p_centr(ceil(length(p_centr)/2))], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii == 3 || ii ==4
        rectangle('position', [off_set(ii)-1, p_centr(1), 2, p_centr(ceil(length(p_centr)/2))- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    end
    
    hold on
    norm_val= max([nanmean(A_hist_pre(this_idx,:)) nanmean(A_hist_post(this_idx,:))]);
    plot((-nanmean(A_hist_pre(this_idx,:))./norm_val)+off_set(ii), p_centr,   'linewidth', 2, 'Color', f_ord(ii,:))
    plot((nanmean(A_hist_post(this_idx,:))./norm_val)+off_set(ii), p_centr,  'linewidth', 2, 'Color', [f_ord(ii,:) .5])
    
    if ii == 1
        %         text(off_set(ii)-.1, 97, '\leftarrow', 'HorizontalAlignment','right', 'Interpreter', 'TeX', 'fontsize', 8)
        text(off_set(ii)-.1, p_centr(end), 'pre', 'HorizontalAlignment','right', 'Interpreter', 'TeX', 'fontsize', 8, 'VerticalAlignment', 'baseline')
        %         text(off_set(ii)+.1, 97, '\rightarrow', 'HorizontalAlignment','left', 'Interpreter', 'TeX', 'fontsize', 8)
        text(off_set(ii)+.1, p_centr(end), 'post', 'HorizontalAlignment','left', 'Interpreter', 'TeX', 'fontsize', 8, 'VerticalAlignment', 'baseline')
        
    end
end
set(gca, 'xtick', off_set, 'ytick', [0 100])
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
ylabel('position (cm)')
axis('square')
title('REM assembly tuning')




subplot(2,4,1)
off_set = 0:2:8;
for ii = 1:5
    if ii == 1
        this_idx = lt1_idx;
    elseif ii == 2
        this_idx = lt5_idx;
    elseif ii == 3
        this_idx = H1_idx;
    elseif ii == 4
        this_idx = H5_idx;
    elseif ii == 5
        this_idx = HS_idx;
    end
    
    if ii == 1 || ii == 2
        rectangle('position', [off_set(ii), p_centr(1), mode(diff(off_set))*.75, p_centr(end)- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii == 5
        rectangle('position', [off_set(ii), p_centr(ceil(length(p_centr)/2)), mode(diff(off_set))*.75, p_centr(end)- p_centr(ceil(length(p_centr)/2))], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii == 3 || ii ==4
        rectangle('position', [off_set(ii), p_centr(1), mode(diff(off_set))*.75, p_centr(ceil(length(p_centr)/2))- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    end
    
    hold on
    norm_val= max([nanmean(A_hist_wake(this_idx,:))]);
    plot((nanmean(A_hist_wake(this_idx,:))./norm_val)+off_set(ii), p_centr,   'linewidth', 2, 'Color', f_ord(ii,:))
    
    
end
set(gca, 'xtick', off_set, 'ytick', [0 100])
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
ylabel('position (cm)')
axis('square')
title('Wake assembly tuning')



%
% figure(99)
% clf
set(gcf, 'Units', 'centimeters', 'Position', [0 0  30 20])
subplot(2,4,2)
% MS_bar_w_err(nWake_A(Cond==1), nWake_A(Cond ==2), c_ord(3,:))
hb = bar([nanmean(nWake_A(Cond==1)), nanmean(nWake_A(Cond ==2)), nanmean(nWake_A(Cond ==3)), nanmean(nWake_A(Cond ==4)), nanmean(nWake_A(Cond ==5))], 'FaceColor', 'flat');
for ii = 1:5
    hb.CData(ii,:) = f_ord(ii,:);
end
hold on
eb = errorbar([nanmean(nWake_A(Cond==1)), nanmean(nWake_A(Cond ==2)), nanmean(nWake_A(Cond ==3)), nanmean(nWake_A(Cond ==4)), nanmean(nWake_A(Cond ==5))], [MS_SEM(nWake_A(Cond==1)), MS_SEM(nWake_A(Cond ==2)), MS_SEM(nWake_A(Cond ==3)), MS_SEM(nWake_A(Cond ==4)), MS_SEM(nWake_A(Cond ==5))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
% ylim(y_lim);
title('Wake Assemblies in Wake')
axis('square')
ylabel('Number of assemblies')


subplot(2,4,6)
hb = bar([nanmean(nPre_A(Cond==1)), nanmean(nPre_A(Cond ==2)), nanmean(nPre_A(Cond ==3)), nanmean(nPre_A(Cond ==4)), nanmean(nPre_A(Cond ==5))] , 'FaceColor', 'flat');
for ii = 1:5
    hb.CData(ii,:) = f_ord(ii,:);
end
hold on

eb = errorbar([nanmean(nPre_A(Cond==1)), nanmean(nPre_A(Cond ==2)), nanmean(nPre_A(Cond ==3)), nanmean(nPre_A(Cond ==4)), nanmean(nPre_A(Cond ==5))], [MS_SEM(nPre_A(Cond==1)), MS_SEM(nPre_A(Cond ==2)), MS_SEM(nPre_A(Cond ==3)), MS_SEM(nPre_A(Cond ==4)), MS_SEM(nPre_A(Cond ==5))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
ylim(y_lim);
title('Pre Assemblies in Wake')
axis('square')


% try to get the pre and post REM assemblies in one plot
wake_mean = [nanmean(nWake_A(Cond==1)), nanmean(nWake_A(Cond ==2)), nanmean(nWake_A(Cond ==3)), nanmean(nWake_A(Cond ==4)), nanmean(nWake_A(Cond ==5))]; 
pre_mean = [nanmean(nPre_sig(Cond==1)), nanmean(nPre_sig(Cond ==2)), nanmean(nPre_sig(Cond ==3)), nanmean(nPre_sig(Cond ==4)), nanmean(nPre_sig(Cond ==5))];
post_mean = [nanmean(nPost_sig(Cond==1)), nanmean(nPost_sig(Cond ==2)), nanmean(nPost_sig(Cond ==3)), nanmean(nPost_sig(Cond ==4)), nanmean(nPost_sig(Cond ==5))];
wake_err = [MS_SEM(nWake_A(Cond==1)), MS_SEM(nWake_A(Cond ==2)), MS_SEM(nWake_A(Cond ==3)), MS_SEM(nWake_A(Cond ==4)), MS_SEM(nWake_A(Cond ==5))];
pre_err = [MS_SEM(nPre_sig(Cond==1)), MS_SEM(nPre_sig(Cond ==2)), MS_SEM(nPre_sig(Cond ==3)), MS_SEM(nPre_sig(Cond ==4)), MS_SEM(nPre_sig(Cond ==5))];
post_err = [MS_SEM(nPost_sig(Cond==1)), MS_SEM(nPost_sig(Cond ==2)), MS_SEM(nPost_sig(Cond ==3)), MS_SEM(nPost_sig(Cond ==4)), MS_SEM(nPost_sig(Cond ==5))];

y = [pre_mean; wake_mean; post_mean]';
err = [pre_err; wake_err; post_err]';

subplot(2,4,7)
cla;
hb = bar(y, 1); % get the bar handles
hold on;
% xpos = ;
pause(.5)

errorbar(hb(1).XData + hb(1).XOffset, y(:,1), err(:,1), 'LineStyle', 'none', ...
    'Color', 'k');
pause(.5)
errorbar(hb(2).XData + hb(2).XOffset, y(:,2), err(:,2), 'LineStyle', 'none', ...
    'Color', 'k');
pause(.5)
errorbar(hb(3).XData + hb(3).XOffset, y(:,3), err(:,3), 'LineStyle', 'none', ...
    'Color', 'k');
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
% ylim(y_lim);
axis('square')
title('Wake assemblies in Pre/Wake/Post REM')
legend({'pre','wake', 'post'}, 'box', 'off')
axis tight
% color

for ii = 1:3
            hb(ii).FaceColor ='flat';
            hb(ii).EdgeColor ='none';

    for jj = 1:5
        hb(ii).CData(jj,:) = f_ord(jj,:);

        % set the alpha for the pre and post. 
        if ii == 1
            hb(ii).FaceAlpha = .7;
        elseif ii == 3
            hb(ii).FaceAlpha = .3;
        end
    end
end

% REACT Str mean(mean Post - mean Pre)
% MS_bar_w_err(nWake_A(Cond==1), nWake_A(Cond ==2), c_ord(3,:))
subplot(2,4,8)
cla
hb = bar([nanmean(ReAct_str(Cond==1)), nanmean(ReAct_str(Cond ==2)), nanmean(ReAct_str(Cond ==3)), nanmean(ReAct_str(Cond ==4)), nanmean(ReAct_str(Cond ==5))], 'FaceColor', 'flat');
for ii = 1:5
    hb.CData(ii,:) = f_ord(ii,:);
end
hold on
eb = errorbar([nanmean(ReAct_str(Cond==1)), nanmean(ReAct_str(Cond ==2)), nanmean(ReAct_str(Cond ==3)), nanmean(ReAct_str(Cond ==4)), nanmean(ReAct_str(Cond ==5))],...
    [MS_SEM(ReAct_str(Cond==1)), MS_SEM(ReAct_str(Cond ==2)), MS_SEM(ReAct_str(Cond ==3)), MS_SEM(ReAct_str(Cond ==4)), MS_SEM(ReAct_str(Cond ==5))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
% ylim(y_lim);
axis('square')
title('Post REM Reactivation Strength')


subplot(2,4,4)
cla
hb = bar([nanmean(nWake_R(Cond==1)), nanmean(nWake_R(Cond ==2)), nanmean(nWake_R(Cond ==3)), nanmean(nWake_R(Cond ==4)), nanmean(nWake_R(Cond ==5))], 'FaceColor', 'flat');
for ii = 1:5
    hb.CData(ii,:) = f_ord(ii,:);
end
hold on
eb = errorbar([nanmean(nWake_R(Cond==1)), nanmean(nWake_R(Cond ==2)), nanmean(nWake_R(Cond ==3)), nanmean(nWake_R(Cond ==4)), nanmean(nWake_R(Cond ==5))],...
    [MS_SEM(nWake_R(Cond==1)), MS_SEM(nWake_R(Cond ==2)), MS_SEM(nWake_R(Cond ==3)), MS_SEM(nWake_R(Cond ==4)), MS_SEM(nWake_R(Cond ==5))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
% ylim(y_lim);
axis('square')
title('wake activation rate (Hz)')


subplot(2,4,3)
cla
hb = bar([nanmean(wake_prct_move(Cond==1)), nanmean(wake_prct_move(Cond ==2)), nanmean(wake_prct_move(Cond ==3)), nanmean(wake_prct_move(Cond ==4)), nanmean(wake_prct_move(Cond ==5))], 'FaceColor', 'flat');
for ii = 1:5
    hb.CData(ii,:) = f_ord(ii,:);
end
hold on
eb = errorbar([nanmean(wake_prct_move(Cond==1)), nanmean(wake_prct_move(Cond ==2)), nanmean(wake_prct_move(Cond ==3)), nanmean(wake_prct_move(Cond ==4)), nanmean(wake_prct_move(Cond ==5))],...
    [MS_SEM(wake_prct_move(Cond==1)), MS_SEM(wake_prct_move(Cond ==2)), MS_SEM(wake_prct_move(Cond ==3)), MS_SEM(wake_prct_move(Cond ==4)), MS_SEM(wake_prct_move(Cond ==5))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
% ylim(y_lim);
axis('square')
title('Percent of time moving')


set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

print(gcf,  [fig_dir filesep   'Stats_summery_SfN_' strrep(num2str(A_out{1}{1}.info.bin), '.', 'p') 's_bin.pdf'], '-dpdf','-r0')

%% only two
figure(88181)
subplot(2,4,7)
% try to get the pre and post REM assemblies in one plot
wake_mean = [nanmean(nWake_A(Cond==1)), nanmean(nWake_A(Cond ==2))]; 
pre_mean = [nanmean(nPre_sig(Cond==1)), nanmean(nPre_sig(Cond ==2))];
post_mean = [nanmean(nPost_sig(Cond==1)), nanmean(nPost_sig(Cond ==2))];
wake_err = [MS_SEM(nWake_A(Cond==1)), MS_SEM(nWake_A(Cond ==2))];
pre_err = [MS_SEM(nPre_sig(Cond==1)), MS_SEM(nPre_sig(Cond ==2))];
post_err = [MS_SEM(nPost_sig(Cond==1)), MS_SEM(nPost_sig(Cond ==2))];

y = [pre_mean;  post_mean]';
err = [pre_err;  post_err]';

cla;
hb = bar(y, 1); % get the bar handles
hold on;
% xpos = ;
pause(.5)

errorbar(hb(1).XData + hb(1).XOffset, y(:,1), err(:,1), 'LineStyle', 'none', ...
    'Color', 'k');
pause(.5)
errorbar(hb(2).XData + hb(2).XOffset, y(:,2), err(:,2), 'LineStyle', 'none', ...
    'Color', 'k');

set(gca, 'XTickLabel', {'Novel', 'Familiar'}, 'XTickLabelRotation', 45)
% ylim(y_lim);
axis('square')
title('Wake assemblies in Pre/Post REM')
legend({'pre' 'post'}, 'box', 'off')
axis tight
% color

for ii = 1:2
    hb(ii).FaceColor ='flat';
    hb(ii).EdgeColor ='none';
    
    for jj = 1:2
        if jj == 1
            hb(ii).CData(jj,:) = f_ord(4,:);
        elseif jj ==2
            hb(ii).CData(jj,:) = f_ord(2,:);
        end
        % set the alpha for the pre and post.
        if ii == 1
            hb(ii).FaceAlpha = .7;
        elseif ii == 2
            hb(ii).FaceAlpha = .3;
        end
    end
end

hold on

set(gca, 'XTickLabel', {'Novel', 'Familiar'}, 'XTickLabelRotation', 45)
ylim([0 25])%max(y, [], 'all')+max(err, [], 'all')*1.5]);
title('wake assemblies in REM')
axis('square')

% Fam vs Anx

subplot(2,4,8)
% try to get the pre and post REM assemblies in one plot
wake_mean = [nanmean(nWake_A(Cond==2)), nanmean(nWake_A(Cond ==3))]; 
pre_mean = [nanmean(nPre_sig(Cond==2)), nanmean(nPre_sig(Cond ==3))];
post_mean = [nanmean(nPost_sig(Cond==2)), nanmean(nPost_sig(Cond ==3))];
wake_err = [MS_SEM(nWake_A(Cond==2)), MS_SEM(nWake_A(Cond ==3))];
pre_err = [MS_SEM(nPre_sig(Cond==2)), MS_SEM(nPre_sig(Cond ==3))];
post_err = [MS_SEM(nPost_sig(Cond==2)), MS_SEM(nPost_sig(Cond ==3))];

y = [pre_mean;  post_mean]';
err = [pre_err;  post_err]';

cla;
hb = bar(y, 1); % get the bar handles
hold on;
% xpos = ;
pause(.5)

errorbar(hb(1).XData + hb(1).XOffset, y(:,1), err(:,1), 'LineStyle', 'none', ...
    'Color', 'k');
pause(.5)
errorbar(hb(2).XData + hb(2).XOffset, y(:,2), err(:,2), 'LineStyle', 'none', ...
    'Color', 'k');

% ylim(y_lim);


for ii = 1:2
    hb(ii).FaceColor ='flat';
    hb(ii).EdgeColor ='none';
    
    for jj = 1:2
        if jj == 1
            hb(ii).CData(jj,:) = f_ord(2,:);
        elseif jj ==2
            hb(ii).CData(jj,:) = f_ord(3,:);
        end
        % set the alpha for the pre and post.
        if ii == 1
            hb(ii).FaceAlpha = .7;
        elseif ii == 2
            hb(ii).FaceAlpha = .3;
        end
    end
end
ylim([0 25])%max(y, [], 'all')+max(err, [], 'all')*1.5]);
set(gca, 'XTickLabel', {'Familiar', 'Anxiety'}, 'XTickLabelRotation', 45)
axis('square')
title('Wake assemblies in REM')
legend({'pre' 'post'}, 'box', 'off')
% color

cond_val = [ones(1,length(nPre_sig(Cond==1))), ones(1,length(nPre_sig(Cond==3)))+1]; 
pre_post = [ones(1,length(cond_val)), ones(1,length(cond_val))+1]; 
cond_val = [cond_val, cond_val]; 

c1 = anovan([nPre_sig(Cond==1)  nPre_sig(Cond ==3), nPost_sig(Cond==1)  nPost_sig(Cond ==3) ], {cond_val, pre_post});
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])



%

subplot(2,4,4)
cla
hb = bar([nanmean(nWake_R(Cond==1)), nanmean(nWake_R(Cond ==2))], 'FaceColor', 'flat');
for ii = 1:2
    hb.CData(ii,:) = f_ord(ii,:);
end
hold on
eb = errorbar([nanmean(nWake_R(Cond==1)), nanmean(nWake_R(Cond ==2))],...
    [MS_SEM(nWake_R(Cond==1)), MS_SEM(nWake_R(Cond ==2))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca, 'XTickLabel', {'Novel', 'Familiar'}, 'XTickLabelRotation', 45)
ylim([0  max([nanmean(nWake_R(Cond==1)), nanmean(nWake_R(Cond ==2))])+max([MS_SEM(nWake_R(Cond==1)), MS_SEM(nWake_R(Cond ==2))])*1.5]);
axis('square')
title('wake activation rate (Hz)')


set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])







print(gcf,  [fig_dir filesep   'Stats_summery_SfN_pre_' strrep(num2str(A_out{1}{1}.info.bin), '.', 'p') 's_bin.pdf'], '-dpdf','-r0')


%% compare the pre and post across novel, familiar and anxiety
figure(19199)
subplot(2,4,1)
% try to get the pre and post REM assemblies in one plot
wake_mean = [nanmean(nWake_A(Cond==1)), nanmean(nWake_A(Cond==2)), nanmean(nWake_A(Cond ==3))]; 
pre_mean = [nanmean(nPre_sig(Cond==1)), nanmean(nPre_sig(Cond==2)), nanmean(nPre_sig(Cond ==3))];
post_mean = [nanmean(nPost_sig(Cond==1)), nanmean(nPost_sig(Cond==2)), nanmean(nPost_sig(Cond ==3))];
wake_err = [MS_SEM(nWake_A(Cond==1)), MS_SEM(nWake_A(Cond==2)), MS_SEM(nWake_A(Cond ==3))];
pre_err = [MS_SEM(nPre_sig(Cond==1)), MS_SEM(nPre_sig(Cond==2)), MS_SEM(nPre_sig(Cond ==3))];
post_err = [MS_SEM(nPost_sig(Cond==1)),MS_SEM(nPost_sig(Cond==2)), MS_SEM(nPost_sig(Cond ==3))];

y = [pre_mean;  post_mean]';
err = [pre_err;  post_err]';

cla;
hb = bar(y, 1); % get the bar handles
hold on;
% xpos = ;
pause(.5)

errorbar(hb(1).XData + hb(1).XOffset, y(:,1), err(:,1), 'LineStyle', 'none', ...
    'Color', 'k');
pause(.5)
errorbar(hb(2).XData + hb(2).XOffset, y(:,2), err(:,2), 'LineStyle', 'none', ...
    'Color', 'k');

% ylim(y_lim);


for ii = 1:2
    hb(ii).FaceColor ='flat';
    hb(ii).EdgeColor ='none';
    
    for jj = 1:3
        if jj == 1
            hb(ii).CData(jj,:) = f_ord(4,:);
        elseif jj ==2
            hb(ii).CData(jj,:) = f_ord(2,:);
        elseif jj ==3
            hb(ii).CData(jj,:) = f_ord(3,:);
        end 
        % set the alpha for the pre and post.
        if ii == 1
            hb(ii).FaceAlpha = .7;
        elseif ii == 2
            hb(ii).FaceAlpha = .3;
        end
    end
end
ylim([0 20])%max(y, [], 'all')+max(err, [], 'all')*1.5]);
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'Anxiety'}, 'XTickLabelRotation', 45)
axis('square')
ylabel({'Num. wake assemblies'; 'reactivated in REM'})
legend({'pre' 'post'}, 'box', 'off')
% color

cond_val = [ones(1,length(nPre_sig(Cond==1))),ones(1,length(nPre_sig(Cond==2)))+1 ones(1,length(nPre_sig(Cond==3)))+2]; 
pre_post = [ones(1,length(cond_val)), ones(1,length(cond_val))+1]; 
cond_val = [cond_val, cond_val]; 

[~,a_tbl, stats] = anovan([nPre_sig(Cond==1)  nPre_sig(Cond ==2) nPre_sig(Cond ==3), nPost_sig(Cond==1)  nPost_sig(Cond ==2) nPost_sig(Cond ==3) ], {cond_val, pre_post}, 2,3,  {'cond', 'pre_post'}, 'off');
a_tbl(:, 4) = []; 

c1 = multcompare(stats,'display','off');
ca_tbl1 = array2table(c1,"VariableNames", ...
    ["A","B","Lower Limit","A-B","Upper Limit","P-value"]);
ca_tbl_cell =["A","B","Lower Limit","A-B","Upper Limit","P-value";  table2cell(ca_tbl1)]; 
ca_tbl_cell(:, 4) = []; 
figure(19199)
subplot(2,4,2)
cla
for ii = 1:size(ca_tbl_cell,1)
    for jj = 1:size(ca_tbl_cell,2)
            if jj == 1 || jj == 2 || jj ==3; 
        f = .75;  
    else
        f = 1;
    end 
        if isnumeric(ca_tbl_cell{ii, jj})
            text((jj-1)*f, -ii-.5, num2str(ca_tbl_cell{ii, jj},3), 'fontsize', 8)
        else
            text((jj-1)*f, -ii, ca_tbl_cell{ii, jj}, 'fontsize',8)
        end
    end
end

for ii = 1:size(a_tbl,1)-2
       for jj = 1:size(a_tbl,2)
        if isnumeric(a_tbl{ii, jj})
            text((jj-1)*f, -ii+size(a_tbl,2)-1.5, num2str(a_tbl{ii, jj},3), 'fontsize', 8)
        else
            text((jj-1)*f, -ii+size(a_tbl,2)-1, a_tbl{ii, jj}, 'fontsize',8)
        end
    end
end
xlim([0 size(ca_tbl_cell, 2)*f])

ylim([-size(ca_tbl_cell,2) size(a_tbl,2)-1])
axis off

%% same thing but for the rate
% try to get the pre and post REM assemblies in one plot
wake_mean = [nanmean(nWake_A(Cond==1)), nanmean(nWake_A(Cond==2)), nanmean(nWake_A(Cond ==3))]; 
pre_mean = [nanmean(nPre_sig(Cond==1)), nanmean(nPre_sig(Cond==2)), nanmean(nPre_sig(Cond ==3))];
post_mean = [nanmean(nPost_sig(Cond==1)), nanmean(nPost_sig(Cond==2)), nanmean(nPost_sig(Cond ==3))];
wake_err = [MS_SEM(nWake_A(Cond==1)), MS_SEM(nWake_A(Cond==2)), MS_SEM(nWake_A(Cond ==3))];
pre_err = [MS_SEM(nPre_sig(Cond==1)), MS_SEM(nPre_sig(Cond==2)), MS_SEM(nPre_sig(Cond ==3))];
post_err = [MS_SEM(nPost_sig(Cond==1)),MS_SEM(nPost_sig(Cond==2)), MS_SEM(nPost_sig(Cond ==3))];


y = [pre_mean; wake_mean; post_mean]';
err = [pre_err; wake_err; post_err]';

y = [pre_mean;  post_mean]';
err = [pre_err;  post_err]';

cla;
hb = bar(y, 1); % get the bar handles
hold on;
% xpos = ;
pause(.5)

errorbar(hb(1).XData + hb(1).XOffset, y(:,1), err(:,1), 'LineStyle', 'none', ...
    'Color', 'k');
pause(.5)
errorbar(hb(2).XData + hb(2).XOffset, y(:,2), err(:,2), 'LineStyle', 'none', ...
    'Color', 'k');

% ylim(y_lim);


for ii = 1:2
    hb(ii).FaceColor ='flat';
    hb(ii).EdgeColor ='none';
    
    for jj = 1:3
        if jj == 1
            hb(ii).CData(jj,:) = f_ord(4,:);
        elseif jj ==2
            hb(ii).CData(jj,:) = f_ord(2,:);
        elseif jj ==3
            hb(ii).CData(jj,:) = f_ord(3,:);
        end 
        % set the alpha for the pre and post.
        if ii == 1
            hb(ii).FaceAlpha = .7;
        elseif ii == 2
            hb(ii).FaceAlpha = .3;
        end
    end
end
ylim([0 20])%max(y, [], 'all')+max(err, [], 'all')*1.5]);
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'Anxiety'}, 'XTickLabelRotation', 45)
axis('square')
ylabel('Wake assemblies in REM')
legend({'pre' 'post'}, 'box', 'off')
% color

cond_val = [ones(1,length(nPre_sig(Cond==1))),ones(1,length(nPre_sig(Cond==2)))+1 ones(1,length(nPre_sig(Cond==3)))+2]; 
pre_post = [ones(1,length(cond_val)), ones(1,length(cond_val))+1]; 
cond_val = [cond_val, cond_val]; 

[~,a_tbl, stats] = anovan([nPre_sig(Cond==1)  nPre_sig(Cond ==2) nPre_sig(Cond ==3), nPost_sig(Cond==1)  nPost_sig(Cond ==2) nPost_sig(Cond ==3) ], {cond_val, pre_post}, 2,3,  {'cond', 'pre_post'}, 'off');
a_tbl(:, 4) = []; 

c1 = multcompare(stats,'display','off');
ca_tbl1 = array2table(c1,"VariableNames", ...
    ["A","B","Lower Limit","A-B","Upper Limit","P-value"]);
ca_tbl_cell =["A","B","Lower Limit","A-B","Upper Limit","P-value";  table2cell(ca_tbl1)]; 
ca_tbl_cell(:, 4) = []; 
figure(19199)
subplot(2,4,2)
cla
for ii = 1:size(ca_tbl_cell,1)
    for jj = 1:size(ca_tbl_cell,2)
            if jj == 1 || jj == 2 || jj ==3; 
        f = .75;  
    else
        f = 1;
    end 
        if isnumeric(ca_tbl_cell{ii, jj})
            text((jj-1)*f, -ii-.5, num2str(ca_tbl_cell{ii, jj},3), 'fontsize', 8)
        else
            text((jj-1)*f, -ii, ca_tbl_cell{ii, jj}, 'fontsize',8)
        end
    end
end

for ii = 1:size(a_tbl,1)-2
       for jj = 1:size(a_tbl,2)
        if isnumeric(a_tbl{ii, jj})
            text((jj-1)*f, -ii+size(a_tbl,2)-1.5, num2str(a_tbl{ii, jj},3), 'fontsize', 8)
        else
            text((jj-1)*f, -ii+size(a_tbl,2)-1, a_tbl{ii, jj}, 'fontsize',8)
        end
    end
end
xlim([0 size(ca_tbl_cell, 2)*f])

ylim([-size(ca_tbl_cell,2) size(a_tbl,2)-1])
axis off

%% plot a specific example
N_idx = []; F_idx = [];
for ii = length(A_out):-1:1
    
    if strcmpi(A_out{ii}{1}.info.subject, 'pv1060') &&  strcmpi(A_out{ii}{1}.info.session, 'HATD1')
        N_idx(ii) = true;
        F_idx(ii) = false;
    elseif strcmpi(A_out{ii}{1}.info.subject, 'pv1060') &&  strcmpi(A_out{ii}{1}.info.session, 'HATD5')
        N_idx(ii) = false;
        F_idx(ii) = true;
    else
        N_idx(ii) = false;
        F_idx(ii) = false;
    end
    
end

F_idx = find(F_idx);
N_idx = find(N_idx);

x_1 = [30 90];
x_2 = [0 60]; 
figure(9908)
clf
for iD = 1:2
    
    if iD == 1
        this_data = A_out{N_idx}{1};
    elseif iD ==2
        this_data = A_out{F_idx}{1};
    end
    
    p_cent = []; p_rate = [];
    
    for ii = length(this_data.P_pos):-1:1
        p_cent(ii) = this_data.map{ii}.cent_z;
        p_rate(ii) = length(this_data.P_loc{ii}.loc_time);
    end
    
    p_rate(p_cent > -1.95) = 0;
    [~, s_idx] = sort(p_rate, 'descend');
    p_rank = s_idx;
    
    if length(this_data.P_pos) >3
        p_idx = p_rank(1:4);
    else
        p_idx = p_rank(1:length(this_data.P_pos));
    end
    
    c_ord = [67, 127, 151; 132 147 35; 255 179 13; 253 22 26]/255;
    
    bx(iD) = subplot(2,2,iD);
    %         yline(this_data.REM_Post_stats.R_thresh, '--', 'color', [.7 .7 .7], 'linewidth', 0.3)
    yline(log10(this_data.REM_Post_stats.R_thresh), '--', 'color', [.7 .7 .7], 'linewidth', 0.3)
    
    hold on
    for ii = length(p_idx):-1:1
        this_proj = this_data.REM_Post_proj(p_idx(ii),:);
        plot(this_data.REM_Post_tvec, log10(this_proj), 'color', c_ord(ii,:), 'linewidth', .5)
        M_max(iD,ii) = max(this_proj);
        
    end
    
    ylim([0 inf])
    y_val = get(gca, 'YTick');
    %         y_val = [0 1 2 2.3979];
    % yTick', [0 1 2 2.3979]
    set(gca, 'YTickLabel', 10.^y_val, 'linewidth', 1);
    
    if iD == 1
        ylabel({'assembly strength'})
    end
    set(gca, 'xtick', []);
    
    if iD == 1
        xlim(x_1);
    elseif iD == 2
        %                 set(gca, 'xtick', 0:60:this_data.REM_Post_tvec(end))
        %         set(gca, 'xticklabel', get(gca, 'xtick') - 90)
        xlim(x_2);
    end
    
    ax2(iD) = subplot(2,2,iD+2);
    for ii = length(p_idx):-1:1
        this_proj_idx = find(this_data.REM_Post_proj(p_idx(ii),:) > this_data.REM_Post_stats.R_thresh);
        line([this_data.REM_Post_tvec(this_proj_idx); this_data.REM_Post_tvec(this_proj_idx)] , [ii-.5; ii+.5], 'color', c_ord(ii,:), 'linewidth', .5)
    end
    
    set(gca, 'ytick', [])
    
    if iD == 1
        xlim(x_1);
        xlabel('time (s)')
        set(gca, 'xtick', [0 60])
        
    elseif iD == 2
        %         set(gca, 'xtick', 0:60:this_data.REM_Post_tvec(end))
        %         set(gca, 'xticklabel', get(gca, 'xtick') - 90)
        xlim(x_2);
        %         x_lim = get(gca, 'xlim');
        set(gca, 'xtick', [])
                                    set(gca, 'xtick', [0 60])
        
    end
    
end
subplot(2,2,1)
linkaxes(bx, 'y')
ylim([0 3])
set(gcf,'PaperUnits','inches', 'Units', 'inches');
set(gcf, 'position', [5 5 7 2])
set(gcf,'PaperSize', [7, 2]);
set(gca,'xlimmode','manual','ylimmode','manual')
print(gcf, '-dpdf', [fig_dir filesep 'Fig3_example.pdf'])

% save an h5 with the data for plotting.
fname = ['SFN_Anx_assembly_.h5'];

if exist(fname, 'file')
    delete(fname)
end

hdf5write(fname, '/mouse', string(A_out{F_idx}{1}.info.subject));
hdf5write(fname, '/N_condition', string(A_out{N_idx}{1}.info.session),'WriteMode', 'append');
hdf5write(fname, '/F_condition', string(A_out{F_idx}{1}.info.session),'WriteMode', 'append');

hdf5write(fname, '/N_REM_post_rThresh', A_out{N_idx}{1}.REM_Post_stats.R_thresh,'WriteMode', 'append');
hdf5write(fname, '/F_REM_post_rThresh', A_out{F_idx}{1}.REM_Post_stats.R_thresh,'WriteMode', 'append');
hdf5write(fname, '/N_REM_post_proj', A_out{N_idx}{1}.REM_Post_proj(p_idx,:),'WriteMode', 'append');
hdf5write(fname, '/F_REM_post_proj', A_out{F_idx}{1}.REM_Post_proj(p_idx,:),'WriteMode', 'append');
hdf5write(fname, '/N_REM_post_tvec', A_out{N_idx}{1}.REM_Post_tvec,'WriteMode', 'append');
hdf5write(fname, '/F_REM_post_tvec', A_out{F_idx}{1}.REM_Post_tvec,'WriteMode', 'append');

hdf5write(fname, '/N_REM_post_xlim', int16([90 150]),'WriteMode', 'append');
hdf5write(fname, '/F_REM_post_xlim', int16([0 60]),'WriteMode', 'append');


%% convert all assemblies into a giant table

A_tbl = []; 
sub_list = []; 
sess_list = [];
win_list = []; 
ReAct_list = []; 
trk_peak = []; 
trk_bias = []; 
maps = []; 

Pre_rate_list = []; 
Post_rate_list = []; 
Pre_str_list = []; 
Post_str_list = []; 


for iA  = 1:length(A_out)
    for iB = 1:length(A_out{iA})
        
        for ii = 1:size(A_out{iA}{iB}.P_pos,2)
            
            sub_list{end+1} =  A_out{iA}{iB}.info.subject; 
            sess_list{end+1} = A_out{iA}{iB}.info.session;
            win_list(end+1) = A_out{iA}{iB}.info.bin; 
%             ReAct_list(end+1) = A_out{iA}{iB}.ReAct(ii); 
            Pre_rate_list(end+1) = A_out{iA}{iB}.REM_Pre_stats.rate(ii);
            Post_rate_list(end+1) = A_out{iA}{iB}.REM_Post_stats.rate(ii);

            
            % compute average expression strength of each assembly. 
            Pre_str_list(end+1) = mean(A_out{iA}{iB}.REM_Pre_proj(ii,:));  
            Post_str_list(end+1) = mean(A_out{iA}{iB}.REM_Post_proj(ii,:));  

            ReAct_list(end+1) = Post_str_list(end) - Pre_str_list(end); 
            % get the peak of the assembly mean map. 
            [~, pk_loc] = findpeaks(A_out{iA}{iB}.Place_map{ii}.map_mean, A_out{iA}{iB}.Place_map{ii}.bins,'SortStr', 'descend');
            pk_loc = pk_loc(1); 
            maps(end+1,:) = A_out{iA}{iB}.Place_map{ii}.map_mean; 
            
            trk_peak(end+1)  = pk_loc; 
            
            % classify the peak into open, closed, or transition on HAT
            % days
            if contains(A_out{iA}{iB}.info.session, {'HATD1', 'HATD3', 'HATD5'})
                
                if pk_loc > 60
                    trk_bias{end+1} = 'O';
                elseif (40 < pk_loc) && (pk_loc< 60)
                    trk_bias{end+1} = 'T';
                elseif pk_loc < 40
                    trk_bias{end+1} = 'C';
                end
          
            elseif contains(A_out{iA}{iB}.info.session, {'HATDSwitch'})
                if pk_loc > 60
                    trk_bias{end+1} = 'C';
                elseif (40 < pk_loc) && (pk_loc< 60)
                    trk_bias{end+1} = 'T';
                elseif pk_loc < 40
                    trk_bias{end+1} = 'O';
                end
            else
                trk_bias{end+1} = 'NaN'; 
            end

        end
        
        
    end
end


A_tbl = table(sub_list',sess_list',win_list',ReAct_list', Pre_rate_list', Post_rate_list',Pre_str_list', Post_str_list', trk_peak',trk_bias',...
    'VariableNames', {'subject','session','Windo_size','ReAct_str','Pre_rate', 'Post_rate','Pre_str', 'Post_str','trk_peak','trk_bias'}); 

A_tbl.session = categorical(A_tbl.session); 
A_tbl.subject = categorical(A_tbl.subject); 
A_tbl.trk_bias = categorical(A_tbl.trk_bias); 

%%


% to do:  scatter of pre vs post str to match van de ven paper. 
figure(9009)
clf
subplot(2,2,1)
hold on
b1 = A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD1'))\A_tbl.Post_rate(ismember(A_tbl.session, 'LTD1')); 
plot(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD1')), A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD1'))*b1,'color', c_ord(1,:), 'linewidth', 3); 
b1 = A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5'))\A_tbl.Post_rate(ismember(A_tbl.session, 'LTD5')); 
plot(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5')), A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5'))*b1, 'color',[.6 .6 .6], 'linewidth', 3);
s1 = scatter(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD1')),A_tbl.Post_rate(ismember(A_tbl.session, 'LTD1')), 20, c_ord(1,:), 'filled');
s2 = scatter(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5')),A_tbl.Post_rate(ismember(A_tbl.session, 'LTD5')), 20,  [.6 .6 .6], 'filled');
 
xlim([0 5])
ylim([0 5])
xl = xline(0);
yline(0)
ylabel('Post rate /min')
xlabel('Pre rate /min')
line([0 5], [0 5], 'color', xl.Color, 'LineWidth', xl.LineWidth, 'linestyle', '--')
legend([s1 s2],'Linear', 'Anxiety', 'box', 'off')

subplot(2,2,2)
hold on
b1 = A_tbl.Pre_str(ismember(A_tbl.session, 'LTD1'))\A_tbl.Post_str(ismember(A_tbl.session, 'LTD1')); 
plot(A_tbl.Pre_str(ismember(A_tbl.session, 'LTD1')), A_tbl.Pre_str(ismember(A_tbl.session, 'LTD1'))*b1,'color', c_ord(1,:), 'linewidth', 3); 
b1 = A_tbl.Pre_str(ismember(A_tbl.session, 'LTD5'))\A_tbl.Post_str(ismember(A_tbl.session, 'LTD5')); 
plot(A_tbl.Pre_str(ismember(A_tbl.session, 'LTD5')), A_tbl.Pre_str(ismember(A_tbl.session, 'LTD5'))*b1, 'color',[.6 .6 .6], 'linewidth', 3); 
s1 = scatter(A_tbl.Pre_str(ismember(A_tbl.session, 'LTD1')),A_tbl.Post_str(ismember(A_tbl.session, 'LTD1')), 20, c_ord(1,:), 'filled');
s2 = scatter(A_tbl.Pre_str(ismember(A_tbl.session, 'LTD5')),A_tbl.Post_str(ismember(A_tbl.session, 'LTD5')), 20, [.6 .6 .6], 'filled');

xlim([-0.2 1.5])
ylim([-.2 1.5])
xl = xline(0);
yline(0)
ylabel('Post str')
xlabel('Pre str')
line([-.2 1.5], [-.2 1.5], 'color', xl.Color, 'LineWidth', xl.LineWidth, 'linestyle', '--')
legend([s1 s2],'Linear', 'Anxiety', 'box', 'off')


subplot(2,2,3)
cla
hold on
b1 = A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5'))\A_tbl.Post_rate(ismember(A_tbl.session, 'LTD5')); 
plot(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5')), A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5'))*b1,'color',  [.6 .6 .6], 'linewidth', 3); 
b1 = A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD1'))\A_tbl.Post_rate(ismember(A_tbl.session, 'HATD1')); 
plot(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD1')), A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD1'))*b1, 'color',c_ord(4,:), 'linewidth', 3); 
s1 = scatter(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5')),A_tbl.Post_rate(ismember(A_tbl.session, 'LTD5')), 20, [.6 .6 .6], 'filled');
s2 = scatter(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD1')),A_tbl.Post_rate(ismember(A_tbl.session, 'HATD1')), 20, c_ord(4,:), 'filled');
xlim([0 5])
ylim([0 5])
xl = xline(0);
yline(0)
ylabel('Post rate /min')
xlabel('Pre rate /min')
line([0 5], [0 5], 'color', xl.Color, 'LineWidth', xl.LineWidth, 'linestyle', '--')
legend([s1 s2],'Linear', 'Anxiety', 'box', 'off')

subplot(2,2,4)
hold on
b1 = A_tbl.Pre_str(ismember(A_tbl.session, 'LTD5'))\A_tbl.Post_str(ismember(A_tbl.session, 'LTD5')); 
plot(A_tbl.Pre_str(ismember(A_tbl.session, 'LTD5')), A_tbl.Pre_str(ismember(A_tbl.session, 'LTD5'))*b1,'color',  [.6 .6 .6], 'linewidth', 3); 
b1 = A_tbl.Pre_str(ismember(A_tbl.session, 'HATD1'))\A_tbl.Post_str(ismember(A_tbl.session, 'HATD1')); 
plot(A_tbl.Pre_str(ismember(A_tbl.session, 'HATD1')), A_tbl.Pre_str(ismember(A_tbl.session, 'HATD1'))*b1, 'color',c_ord(4,:), 'linewidth', 3); 
s1 = scatter(A_tbl.Pre_str(ismember(A_tbl.session, 'LTD5')),A_tbl.Post_str(ismember(A_tbl.session, 'LTD5')), 20, [.6 .6 .6], 'filled');
s2 = scatter(A_tbl.Pre_str(ismember(A_tbl.session, 'HATD1')),A_tbl.Post_str(ismember(A_tbl.session, 'HATD1')), 20, c_ord(4,:), 'filled');


xlim([-0.2 1.5])
ylim([-.2 1.5])
xl = xline(0); 
yline(0)
line([-.2 1.5], [-.2 1.5], 'color', xl.Color, 'LineWidth', xl.LineWidth, 'linestyle', '--')
ylabel('Post str')
xlabel('Pre str')
legend([s1 s2],'Linear', 'Anxiety', 'box', 'off')

%%
figure(99999)
clf
subplot(2,3,1)
hold on
% open_idx = (ismember(A_tbl.trk_bias, 'O')); 
% scatter(A_tbl.session(open_idx), A_tbl.ReAct_str(open_idx), 55, A_tbl.subject(open_idx), 'o')
% 
% close_idx = (ismember(A_tbl.trk_bias, 'C')); 
% scatter(A_tbl.session(close_idx), A_tbl.ReAct_str(close_idx), 55, A_tbl.subject(close_idx), 'x')
% 
% trans_idx = (ismember(A_tbl.trk_bias, 'T')); 
% scatter(A_tbl.session(trans_idx), A_tbl.ReAct_str(trans_idx), 55, A_tbl.subject(trans_idx), 's')
% 
% non_idx = ~open_idx | ~close_idx | ~trans_idx; 
% 
% scatter(A_tbl.session(non_idx), A_tbl.ReAct_str(non_idx), 55, A_tbl.subject(non_idx), '.')

    %%%%%%% pre rate

     sc{1} = scatter(1+sort(MS_randn_range(length(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD1'))), 1, -.2, .2)), A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD1')),25, c_ord(1,:), 'filled');
               line([.9 1.1], [mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD1'))) mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD1')))], 'color',  'k', 'linewidth', 2)
     eb = errorbar(1, mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD1'))), MS_SEM(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD1'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
     sc{1} = scatter(2+sort(MS_randn_range(length(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD5'))), 1, -.1, .1)), A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD5')),25, [.8 .8 .8], 'filled');
            line([1.9 2.1], [mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD5'))) mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD5')))], 'color',  'k', 'linewidth', 2)
     eb = errorbar(2, mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD5'))), MS_SEM(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD5'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
    sc{1} = scatter(3+sort(MS_randn_range(length(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1'))), 1, -.1, .1)), A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1')),25, c_ord(4,:), 'filled');
            line([2.9 3.1], [mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1'))) mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1')))], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3, mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1'))), MS_SEM(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
    sc{1} = scatter(4+sort(MS_randn_range(length(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD5'))), 1, -.1, .1)), A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD5')),25, c_ord(3,:), 'filled');
        line([3.9 4.1], [mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD5'))) mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD5')))], 'color',  'k', 'linewidth', 2)
    eb = errorbar(4, mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD5'))), MS_SEM(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD5'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
    sc{1} = scatter(5+sort(MS_randn_range(length(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATDSwitch'))), 1, -.1, .1)), A_tbl.ReAct_str(ismember(A_tbl.session, 'HATDSwitch')),25, c_ord(2,:), 'filled');
    line([4.9 5.1], [mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATDSwitch'))) mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATDSwitch')))], 'color',  'k', 'linewidth', 2)
    eb = errorbar(5, mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATDSwitch'))), MS_SEM(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATDSwitch'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
%         ylabel({'Pre REM assembly rate' ; 'activations/mins'})

yline(0)
ylabel({'Reactivation strength'; 'weaker <-    |    -> stronger'})
    set(gca, 'xtick', 1:5, 'xticklabels', {'LTD1', 'LTD5', 'HATD1', 'HATD5', 'HATDS'})
    ylim([-1 1])

% same plot but easier to read
subplot(2,3,4)
cla

hold on

% bar(1:2, [mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD1'))); mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD2')))])
% HAT1_ReAct = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1')); 

% x_data = 
ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD1'));
    sc{1} = scatter(1+sort(MS_randn_range(length(ydata), 1, -.1, .1)), ydata,25, c_ord(1,:), 'filled');
    line([.9 1.1], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(1, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    
ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD5'));
    sc{1} = scatter(2+sort(MS_randn_range(length(ydata), 1, -.1, .1)), ydata,25, [.8 .8 .8], 'filled');
    line([1.9 2.1], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(2, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
   
    % hat varients
ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1') & ismember(A_tbl.trk_bias, 'O'));
    sc{1} = scatter(3+sort(MS_randn_range(length(ydata), 1, -.2, -.1)), ydata,25, c_ord(4,:));
    line([2.8 2.9], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(2.85, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1') & ismember(A_tbl.trk_bias, 'T'));
    sc{1} = scatter(3+sort(MS_randn_range(length(ydata), 1, -0.05, .05)), ydata,25, c_ord(4,:)*.5,'x');
    line([2.95 3.05], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1') & ismember(A_tbl.trk_bias, 'C'));
    sc{1} = scatter(3+sort(MS_randn_range(length(ydata), 1, .1, .2)), ydata,25, c_ord(4,:)*.8,'o', 'filled');
    line([3.1 3.2], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3.15, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    %hat 5
ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD5') & ismember(A_tbl.trk_bias, 'O'));
    sc{1} = scatter(4+sort(MS_randn_range(length(ydata), 1, -.2, -.1)), ydata,25, c_ord(3,:));
    line([3.8 3.9], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3.85, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD5') & ismember(A_tbl.trk_bias, 'T'));
    sc{1} = scatter(4+sort(MS_randn_range(length(ydata), 1, -0.05, .05)), ydata,25, c_ord(3,:)*.5,'x');
    line([3.95 4.05], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(4, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD5') & ismember(A_tbl.trk_bias, 'C'));
    sc{1} = scatter(4+sort(MS_randn_range(length(ydata), 1, .1, .2)), ydata,25, c_ord(3,:)*.8,'o', 'filled');
    line([4.1 4.2], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(4.15, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    
        %hat switch

ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATDSwitch') & ismember(A_tbl.trk_bias, 'O'));
    sc{1} = scatter(5+sort(MS_randn_range(length(ydata), 1, -.2, -.1)), ydata,25, c_ord(2,:));
    line([4.8 4.9], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(4.85, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATDSwitch') & ismember(A_tbl.trk_bias, 'T'));
    sc{1} = scatter(5+sort(MS_randn_range(length(ydata), 1, -0.05, .05)), ydata,25, c_ord(2,:)*.5,'x');
    line([4.95 5.05], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(5, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATDSwitch') & ismember(A_tbl.trk_bias, 'C'));
    sc{1} = scatter(5+sort(MS_randn_range(length(ydata), 1, .1, .2)), ydata,25, c_ord(2,:)*.8,'o', 'filled');
    line([5.1 5.2], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(5.15, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    
    yline(0)
ylabel({'Reactivation strength'; 'weaker <-    |    -> stronger'})

       set(gca, 'xtick', 1:5, 'xticklabels', {'LTD1', 'LTD5', 'HATD1', 'HATD5', 'HATDS'})
        ylim([-1 1])

%%%%%%% pre rate
    %%%%%%% pre rate
    subplot(2,3,2)
    cla
    hold on
     sc{1} = scatter(1+sort(MS_randn_range(length(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD1'))), 1, -.2, .2)), A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD1')),25, c_ord(1,:), 'filled');
               line([.9 1.1], [mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD1'))) mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD1')))], 'color',  'k', 'linewidth', 2)
     eb = errorbar(1, mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD1'))), MS_SEM(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD1'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
     sc{1} = scatter(2+sort(MS_randn_range(length(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5'))), 1, -.1, .1)), A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5')),25, [.8 .8 .8], 'filled');
            line([1.9 2.1], [mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5'))) mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5')))], 'color',  'k', 'linewidth', 2)
     eb = errorbar(2, mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5'))), MS_SEM(A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
    sc{1} = scatter(3+sort(MS_randn_range(length(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD1'))), 1, -.1, .1)), A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD1')),25, c_ord(4,:), 'filled');
            line([2.9 3.1], [mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD1'))) mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD1')))], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3, mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD1'))), MS_SEM(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD1'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
    sc{1} = scatter(4+sort(MS_randn_range(length(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD5'))), 1, -.1, .1)), A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD5')),25, c_ord(3,:), 'filled');
        line([3.9 4.1], [mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD5'))) mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD5')))], 'color',  'k', 'linewidth', 2)
    eb = errorbar(4, mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD5'))), MS_SEM(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD5'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
    sc{1} = scatter(5+sort(MS_randn_range(length(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATDSwitch'))), 1, -.1, .1)), A_tbl.Pre_rate(ismember(A_tbl.session, 'HATDSwitch')),25, c_ord(2,:), 'filled');
    line([4.9 5.1], [mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATDSwitch'))) mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATDSwitch')))], 'color',  'k', 'linewidth', 2)
    eb = errorbar(5, mean(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATDSwitch'))), MS_SEM(A_tbl.Pre_rate(ismember(A_tbl.session, 'HATDSwitch'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
        ylabel({'Pre REM assembly rate' ; 'activations/mins'})
    set(gca, 'xtick', 1:5, 'xticklabels', {'LTD1', 'LTD5', 'HATD1', 'HATD5', 'HATDS'})

subplot(2,3,5)
cla

hold on

% bar(1:2, [mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD1'))); mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD2')))])
% HAT1_ReAct = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1')); 

% x_data = 
ydata = A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD1'));
    sc{1} = scatter(1+sort(MS_randn_range(length(ydata), 1, -.1, .1)), ydata,25, c_ord(1,:), 'filled');
    line([.9 1.1], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(1, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    
ydata = A_tbl.Pre_rate(ismember(A_tbl.session, 'LTD5'));
    sc{1} = scatter(2+sort(MS_randn_range(length(ydata), 1, -.1, .1)), ydata,25, [.8 .8 .8], 'filled');
    line([1.9 2.1], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(2, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
   
    % hat varients
ydata = A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD1') & ismember(A_tbl.trk_bias, 'O'));
    sc{1} = scatter(3+sort(MS_randn_range(length(ydata), 1, -.2, -.1)), ydata,25, c_ord(4,:));
    line([2.8 2.9], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(2.85, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    ydata = A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD1') & ismember(A_tbl.trk_bias, 'T'));
    sc{1} = scatter(3+sort(MS_randn_range(length(ydata), 1, -0.05, .05)), ydata,25, c_ord(4,:)*.5,'x');
    line([2.95 3.05], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD1') & ismember(A_tbl.trk_bias, 'C'));
    sc{1} = scatter(3+sort(MS_randn_range(length(ydata), 1, .1, .2)), ydata,25, c_ord(4,:)*.8,'o', 'filled');
    line([3.1 3.2], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3.15, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    %hat 5
ydata = A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD5') & ismember(A_tbl.trk_bias, 'O'));
    sc{1} = scatter(4+sort(MS_randn_range(length(ydata), 1, -.2, -.1)), ydata,25, c_ord(3,:));
    line([3.8 3.9], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3.85, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD5') & ismember(A_tbl.trk_bias, 'T'));
    sc{1} = scatter(4+sort(MS_randn_range(length(ydata), 1, -0.05, .05)), ydata,25, c_ord(3,:)*.5,'x');
    line([3.95 4.05], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(4, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.Pre_rate(ismember(A_tbl.session, 'HATD5') & ismember(A_tbl.trk_bias, 'C'));
    sc{1} = scatter(4+sort(MS_randn_range(length(ydata), 1, .1, .2)), ydata,25, c_ord(3,:)*.8,'o', 'filled');
    line([4.1 4.2], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(4.15, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    
        %hat switch
ydata = A_tbl.Pre_rate(ismember(A_tbl.session, 'HATDSwitch') & ismember(A_tbl.trk_bias, 'O'));
    sc{1} = scatter(5+sort(MS_randn_range(length(ydata), 1, -.2, -.1)), ydata,25, c_ord(2,:));
    line([4.8 4.9], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(4.85, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.Pre_rate(ismember(A_tbl.session, 'HATDSwitch') & ismember(A_tbl.trk_bias, 'T'));
    sc{1} = scatter(5+sort(MS_randn_range(length(ydata), 1, -0.05, .05)), ydata,25, c_ord(2,:)*.5,'x');
    line([4.95 5.05], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(5, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.Pre_rate(ismember(A_tbl.session, 'HATDSwitch') & ismember(A_tbl.trk_bias, 'C'));
    sc{1} = scatter(5+sort(MS_randn_range(length(ydata), 1, .1, .2)), ydata,25, c_ord(2,:)*.8,'o', 'filled');
    line([5.1 5.2], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(5.15, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    
    ylabel({'Pre REM assembly rate' ; 'activations/mins'})
    set(gca, 'xtick', 1:5, 'xticklabels', {'LTD1', 'LTD5', 'HATD1', 'HATD5', 'HATDS'})
    
    
    
    %%%%%%% post rate
    subplot(2,3,3)
    cla
    hold on
     sc{1} = scatter(1+sort(MS_randn_range(length(A_tbl.Post_rate(ismember(A_tbl.session, 'LTD1'))), 1, -.2, .2)), A_tbl.Post_rate(ismember(A_tbl.session, 'LTD1')),25, c_ord(1,:), 'filled');
               line([.9 1.1], [mean(A_tbl.Post_rate(ismember(A_tbl.session, 'LTD1'))) mean(A_tbl.Post_rate(ismember(A_tbl.session, 'LTD1')))], 'color',  'k', 'linewidth', 2)
     eb = errorbar(1, mean(A_tbl.Post_rate(ismember(A_tbl.session, 'LTD1'))), MS_SEM(A_tbl.Post_rate(ismember(A_tbl.session, 'LTD1'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
     sc{1} = scatter(2+sort(MS_randn_range(length(A_tbl.Post_rate(ismember(A_tbl.session, 'LTD5'))), 1, -.1, .1)), A_tbl.Post_rate(ismember(A_tbl.session, 'LTD5')),25, [.8 .8 .8], 'filled');
            line([1.9 2.1], [mean(A_tbl.Post_rate(ismember(A_tbl.session, 'LTD5'))) mean(A_tbl.Post_rate(ismember(A_tbl.session, 'LTD5')))], 'color',  'k', 'linewidth', 2)
     eb = errorbar(2, mean(A_tbl.Post_rate(ismember(A_tbl.session, 'LTD5'))), MS_SEM(A_tbl.Post_rate(ismember(A_tbl.session, 'LTD5'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
    sc{1} = scatter(3+sort(MS_randn_range(length(A_tbl.Post_rate(ismember(A_tbl.session, 'HATD1'))), 1, -.1, .1)), A_tbl.Post_rate(ismember(A_tbl.session, 'HATD1')),25, c_ord(4,:), 'filled');
            line([2.9 3.1], [mean(A_tbl.Post_rate(ismember(A_tbl.session, 'HATD1'))) mean(A_tbl.Post_rate(ismember(A_tbl.session, 'HATD1')))], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3, mean(A_tbl.Post_rate(ismember(A_tbl.session, 'HATD1'))), MS_SEM(A_tbl.Post_rate(ismember(A_tbl.session, 'HATD1'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
    sc{1} = scatter(4+sort(MS_randn_range(length(A_tbl.Post_rate(ismember(A_tbl.session, 'HATD5'))), 1, -.1, .1)), A_tbl.Post_rate(ismember(A_tbl.session, 'HATD5')),25, c_ord(3,:), 'filled');
        line([3.9 4.1], [mean(A_tbl.Post_rate(ismember(A_tbl.session, 'HATD5'))) mean(A_tbl.Post_rate(ismember(A_tbl.session, 'HATD5')))], 'color',  'k', 'linewidth', 2)
    eb = errorbar(4, mean(A_tbl.Post_rate(ismember(A_tbl.session, 'HATD5'))), MS_SEM(A_tbl.Post_rate(ismember(A_tbl.session, 'HATD5'))), 'color', 'k'); 
    eb.LineWidth = 2;
    
    sc{1} = scatter(5+sort(MS_randn_range(length(A_tbl.Post_rate(ismember(A_tbl.session, 'HATDSwitch'))), 1, -.1, .1)), A_tbl.Post_rate(ismember(A_tbl.session, 'HATDSwitch')),25, c_ord(2,:), 'filled');
    line([4.9 5.1], [mean(A_tbl.Post_rate(ismember(A_tbl.session, 'HATDSwitch'))) mean(A_tbl.Post_rate(ismember(A_tbl.session, 'HATDSwitch')))], 'color',  'k', 'linewidth', 2)
    eb = errorbar(5, mean(A_tbl.Post_rate(ismember(A_tbl.session, 'HATDSwitch'))), MS_SEM(A_tbl.Post_rate(ismember(A_tbl.session, 'HATDSwitch'))), 'color', 'k'); 
    eb.LineWidth = 2;

        ylabel({'Post REM assembly rate' ; 'activations/mins'})
    set(gca, 'xtick', 1:5, 'xticklabels', {'LTD1', 'LTD5', 'HATD1', 'HATD5', 'HATDS'})
    
subplot(2,3,6)
cla

hold on

% bar(1:2, [mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD1'))); mean(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD2')))])
% HAT1_ReAct = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1')); 

% x_data = 
ydata = A_tbl.Post_rate(ismember(A_tbl.session, 'LTD1'));
    sc{1} = scatter(1+sort(MS_randn_range(length(ydata), 1, -.1, .1)), ydata,25, c_ord(1,:), 'filled');
    line([.9 1.1], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(1, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    
ydata = A_tbl.Post_rate(ismember(A_tbl.session, 'LTD5'));
    sc{1} = scatter(2+sort(MS_randn_range(length(ydata), 1, -.1, .1)), ydata,25, [.8 .8 .8], 'filled');
    line([1.9 2.1], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(2, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
   
    % hat varients
ydata = A_tbl.Post_rate(ismember(A_tbl.session, 'HATD1') & ismember(A_tbl.trk_bias, 'O'));
    sc{1} = scatter(3+sort(MS_randn_range(length(ydata), 1, -.2, -.1)), ydata,25, c_ord(4,:));
    line([2.8 2.9], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(2.85, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    ydata = A_tbl.Post_rate(ismember(A_tbl.session, 'HATD1') & ismember(A_tbl.trk_bias, 'T'));
    sc{1} = scatter(3+sort(MS_randn_range(length(ydata), 1, -0.05, .05)), ydata,25, c_ord(4,:)*.5,'x');
    line([2.95 3.05], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.Post_rate(ismember(A_tbl.session, 'HATD1') & ismember(A_tbl.trk_bias, 'C'));
    sc{1} = scatter(3+sort(MS_randn_range(length(ydata), 1, .1, .2)), ydata,25, c_ord(4,:)*.8,'o', 'filled');
    line([3.1 3.2], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3.15, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    %hat 5
ydata = A_tbl.Post_rate(ismember(A_tbl.session, 'HATD5') & ismember(A_tbl.trk_bias, 'O'));
    sc{1} = scatter(4+sort(MS_randn_range(length(ydata), 1, -.2, -.1)), ydata,25, c_ord(3,:));
    line([3.8 3.9], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3.85, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.Post_rate(ismember(A_tbl.session, 'HATD5') & ismember(A_tbl.trk_bias, 'T'));
    sc{1} = scatter(4+sort(MS_randn_range(length(ydata), 1, -0.05, .05)), ydata,25, c_ord(3,:)*.5,'x');
    line([3.95 4.05], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(4, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.Post_rate(ismember(A_tbl.session, 'HATD5') & ismember(A_tbl.trk_bias, 'C'));
    sc{1} = scatter(4+sort(MS_randn_range(length(ydata), 1, .1, .2)), ydata,25, c_ord(3,:)*.8,'o', 'filled');
    line([4.1 4.2], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(4.15, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    
        %hat switch
ydata = A_tbl.Post_rate(ismember(A_tbl.session, 'HATDSwitch') & ismember(A_tbl.trk_bias, 'O'));
    sc{1} = scatter(5+sort(MS_randn_range(length(ydata), 1, -.2, -.1)), ydata,25, c_ord(2,:));
    line([4.8 4.9], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(4.85, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.Post_rate(ismember(A_tbl.session, 'HATDSwitch') & ismember(A_tbl.trk_bias, 'T'));
    sc{1} = scatter(5+sort(MS_randn_range(length(ydata), 1, -0.05, .05)), ydata,25, c_ord(2,:)*.5,'x');
    line([4.95 5.05], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(5, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.Post_rate(ismember(A_tbl.session, 'HATDSwitch') & ismember(A_tbl.trk_bias, 'C'));
    sc{1} = scatter(5+sort(MS_randn_range(length(ydata), 1, .1, .2)), ydata,25, c_ord(2,:)*.8,'o', 'filled');
    line([5.1 5.2], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(5.15, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    
    ylabel({'Post REM assembly rate' ; 'activations/mins'})
    set(gca, 'xtick', 1:5, 'xticklabels', {'LTD1', 'LTD5', 'HATD1', 'HATD5', 'HATDS'})
%     sc{1} = scatter(2+sort(MS_randn_range(length(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD5'))), 1, -.1, .1)), A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD5')),25, c_ord(1,:), 'filled');
%     sc{1} = scatter(3+sort(MS_randn_range(length(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1'))), 1, -.1, .1)), A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1')),25, c_ord(1,:), 'filled');
%     sc{1} = scatter(4+sort(MS_randn_range(length(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD5'))), 1, -.1, .1)), A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD5')),25, c_ord(1,:), 'filled');
%     sc{1} = scatter(5+sort(MS_randn_range(length(A_tbl.ReAct_str(ismember(A_tbl.session, 'HATDSwitch'))), 1, -.1, .1)), A_tbl.ReAct_str(ismember(A_tbl.session, 'HATDSwitch')),25, c_ord(1,:), 'filled');





%% stats for population measures
disp('ReAct Str')
[p, tbl, stats] = anovan(A_tbl.ReAct_str, {A_tbl.session}, 'varnames', 'session'); 
c1 = multcompare(stats);
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


% pre & post rate
disp('Pre rate')
[p, tbl, stats] = anovan(A_tbl.Post_rate, {A_tbl.session}, 'varnames', 'session'); 
c1 = multcompare(stats);
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


disp('Post rate')
[p, tbl, stats] = anovan(A_tbl.Post_rate, {A_tbl.session}, 'varnames', 'session'); 
c1 = multcompare(stats);
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

%% 
%% 2 way only/ 

keep_idx = ismember(A_tbl.session, {'LTD5', 'HATD1'});
Rate_lg = [A_tbl.Pre_rate(keep_idx); A_tbl.Post_rate(keep_idx)]; 

Sess_lg = categorical([A_tbl.session(keep_idx); A_tbl.session(keep_idx)]'); 
Sess_lg = removecats(Sess_lg, {'HATDSwitch', 'HATD5', 'LTD1'});
Sess_lg = reordercats(Sess_lg, flipud(categories(Sess_lg)));

Pre_post_lg = categorical([repmat({'Pre'}, 1, length(A_tbl.Pre_rate(keep_idx))) repmat({'Post'}, 1, length(A_tbl.Pre_rate(keep_idx)))])'; %logical([ones(1,length(A_tbl.Pre_rate)), ones(1,length(A_tbl.Post_rate))+1])'; 
Pre_post_lg = reordercats(Pre_post_lg, flipud(categories(Pre_post_lg)));

% keep_idx = ismember(Sess_lg, {'LTD5', 'HATD1'}); 
% can deal with unbalanced. 
[p, tbl, stats] = anovan(Rate_lg, {Sess_lg, Pre_post_lg},'model', 2, 'varnames', {'session', 'Pre-Post'}); % , 
% 
% % [p, tbl, stats] = anovan(A_tbl.ReAct_str, {A_tbl.session,pre_post}, 2,2, {'cond', 'pre_post'});
% % 
c1 = multcompare(stats);
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])

figure(1919)
clf
subplot(2,1,1)
bp = boxchart(Sess_lg, Rate_lg, 'GroupbyColor', Sess_lg);
bp(1).BoxFaceColor = [.6 .6 .6];
bp(2).BoxFaceColor = c_ord(4,:);
ylabel({'REM assembly' ; 'activation rate /min'})
% legend

subplot(2,1,2)
bp = boxchart(Sess_lg, Rate_lg,  'GroupbyColor', Pre_post_lg);
bp(1).BoxFaceColor = [.6 .6 .6];
bp(2).BoxFaceColor = c_ord(4,:);
ylabel({'REM assembly' ; 'activation rate /min'})

legend
%% try repeated measures using just anova (assumes balance)
% A_lg = table(Rate_lg, Sess_lg, Pre_post_lg, 'VariableNames', {'Rate', 'Sess', 'Pre_post'}); 
% 
LA_tbl = A_tbl; 
LA_tbl(~ismember(LA_tbl.session, {'HATD1', 'LTD5'}),:) = []; 

LA_tbl.session = removecats(LA_tbl.session, {'LTD1', 'HATD5', 'HATDSwitch'}); 
% 
% rate_rm = fitrm(LA_tbl, 'Pre_rate - Post_rate ~ session ');
% figure(101)
% clf
% plot(rate_rm)
% ranova_tbl = ranova(rate_rm, 'WithinModel', 'Pre_post');
% % results
% disp(ranova_tbl);
%     multcompare(ranova_tbl)
    % same but with LMM

    
    lg_tbl = table(Rate_lg, Sess_lg', Pre_post_lg, 'VariableNames', {'Rate','Sess', 'Pre_post'}); 
    lg_tbl.Sess = reordercats(lg_tbl.Sess, flipud(categories(lg_tbl.Sess))); 
    lg_tbl.Pre_post = reordercats(lg_tbl.Pre_post, flipud(categories(lg_tbl.Pre_post))); 

%     lg_tbl = sortrows(lg_tbl,{'Sess','Pre_post'},{'descend','descend'});
    
    figure(9898); clf
boxchart(lg_tbl.Sess, lg_tbl.Rate, 'GroupbyColor', lg_tbl.Pre_post)
legend
% set(gca, 'xticklabel', {'Familiar', 'Anxiety'})

rate_lmm = fitlm(lg_tbl, 'Rate ~ Sess * Pre_post')
plot(rate_lmm)
    
    %% try a different way.
    
    
    
%%
figure(8810)
clf
subplot(2,2,2)
boxchart(Rate_lg(keep_idx), A_tbl.Pre_str)


% model
day_mdl = fitrm(A_tbl, 'ReAct_str~sesson')%, 'WithinDesign',table([1 2 3 4 5 6],'VariableNames', "session"));
% repeated measures ANOVA
ranova_tbl = ranova(day_mdl, 'WithinModel', 'Days');
% results
disp(ranova_tbl);
    multcompare(ranova_tbl)




%% simple Linear vs Anixety plots


figure(10000)
clf

subplot(1,3,1)

 [hb, eb, sc, p, stats] = MS_bar_w_err(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD5'))', A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD5'))', [[.6 .6 .6]; c_ord(4,:)], 2, 'ttest2');
 fprintf('LTD5 vs HATD1 ReAct_str ttest: df: %2.0f  tstat: %2.2f p: %2.2f\n', stats.df, stats.tstat, p)
xlim([.5 2.5]);
ylabel({'Post REM reactivation strength'; 'weaker <-    |    -> stronger'})
set(gca, 'xtick', 1:2, 'xticklabel', {'Linear', 'Anxiety'})
yline(0)
 [hb, eb, sc, p, stats] = MS_bar_w_err(A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD5'))', A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD5'))', [[.6 .6 .6]; c_ord(4,:)], 2, 'ttest2');
 fprintf('LTD5 vs HATD1 ReAct_str ttest: df: %2.0f  tstat: %2.2f p: %2.2f\n', stats.df, stats.tstat, p)
xlim([.5 2.5]);
ylabel({'Post REM reactivation strength'; 'weaker <-    |    -> stronger'})
set(gca, 'xtick', 1:2, 'xticklabel', {'Linear', 'Anxiety'})


subplot(1,3,2)

    hold on
% ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'LTD5'));
%     sc{1} = scatter(2+sort(MS_randn_range(length(ydata), 1, -.1, .1)), ydata,25, [.8 .8 .8], 'filled');
%     line([1.9 2.1], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
%     eb = errorbar(2, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
%     eb.LineWidth = 2; 
   
    % hat varients
ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1') & ismember(A_tbl.trk_bias, 'O'));
    sc{1} = scatter(3+sort(MS_randn_range(length(ydata), 1, -.2, -.1)), ydata,25, c_ord(4,:));
    line([2.8 2.9], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(2.85, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
    ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1') & ismember(A_tbl.trk_bias, 'T'));
    sc{1} = scatter(3+sort(MS_randn_range(length(ydata), 1, -0.05, .05)), ydata,25, c_ord(4,:)*.5,'x');
    line([2.95 3.05], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ydata = A_tbl.ReAct_str(ismember(A_tbl.session, 'HATD1') & ismember(A_tbl.trk_bias, 'C'));
    sc{1} = scatter(3+sort(MS_randn_range(length(ydata), 1, .1, .2)), ydata,25, c_ord(4,:)*.8,'o', 'filled');
    line([3.1 3.2], [mean(ydata) mean(ydata)], 'color',  'k', 'linewidth', 2)
    eb = errorbar(3.15, mean(ydata), MS_SEM(ydata), 'color', 'k'); 
    eb.LineWidth = 2; 
    
ylabel({'Post REM reactivation strength'; 'weaker <-    |    -> stronger'})
set(gca, 'xtick', 2:3, 'xticklabel', {'Linear', 'Anxiety'})
xlim([2.75 3.25])
ylim([-.6 .6])
yline(0); 

% one way anova for differences in HATD1 trk bias
disp('HAT Bias ReAct_str')
t_tbl = A_tbl; 
t_tbl.trk_bias = removecats(t_tbl.trk_bias, {'NaN'}); 
[p, tbl, stats] = anovan(t_tbl.ReAct_str(ismember(t_tbl.session, 'HATD1')), {t_tbl.trk_bias(ismember(t_tbl.session, 'HATD1'))}, 'varnames', 'trk_bias'); 

 fprintf('HATD1 trk bias ReAct_str %s: F(%.0f): %2.2f p= %2.2f\n',stats.source, tbl{2,3}, tbl{2,6}, p)
 
c1 = multcompare(stats,'display','off');
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


subplot(1,3,3)
bp = boxchart(Sess_lg, Rate_lg,  'GroupbyColor', Pre_post_lg,'Notch','on');
bp(1).BoxFaceColor = [.6 .6 .6];
bp(2).BoxFaceColor = c_ord(4,:);
ylabel({'REM assembly' ; 'activation rate /min'})

[p, tbl, stats] = anovan(Rate_lg, {Sess_lg, Pre_post_lg},'model', 2, 'varnames', {'session', 'Pre-Post'}); % ,
for ii = 1:length(stats.terms)
 fprintf('HATD1 trk bias ReAct_str %s: F(%.0f): %2.2f p= %2.2f\n',tbl{ii+1,1}, tbl{ii+1,3}, tbl{ii+1,6}, tbl{ii+1,7})
end

% % [p, tbl, stats] = anovan(A_tbl.ReAct_str, {A_tbl.session,pre_post}, 2,2, {'cond', 'pre_post'});
% % 
c1 = multcompare(stats,'display','off');
tbl1 = array2table(c1,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"])


%% wake normalized tuning


% f_ord = [[253, 22, 26];[255, 179, 13]; [132, 147, 35];; [0, 42, 95]]/255;
f_ord = [0.2627 0.4980 0.5922; .7 .7 .7; 0.9922 0.0863 0.1020;1.0000 0.7020 0.0510; 0.5176 0.5765 0.1373]; 
f_post = [f_ord ,[.5, .5, .5, .5, .5]'];
off_set = 0:2:8;

figure(8119)
clf
set(gcf, 'Units', 'centimeters', 'Position', [0 0  30 20])

% hold the histograms for normalization

n_pre_hist = A_hist_pre; 
n_post_hist = A_hist_post; 
n_wake_hist = A_hist_wake; 
n_wake_hist(n_wake_hist == 0) = eps; 
n_pre_hist = n_pre_hist./n_wake_hist;
n_post_hist = n_post_hist./n_wake_hist;

% plot the pre and post REM assembly centroid reactivations.
subplot(1,2,1)
off_set = 0:2.5:10;
for ii = 1:5
    if ii == 1
        this_idx = lt1_idx;
    elseif ii == 2
        this_idx = lt5_idx;
    elseif ii == 3
        this_idx = H1_idx;
    elseif ii == 4
        this_idx = H5_idx;
    elseif ii == 5
        this_idx = HS_idx;
    end
    % add 'walls'
    if ii == 1 || ii == 2
        rectangle('position', [off_set(ii)-1, p_centr(1), 2, p_centr(end)- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii == 5
        rectangle('position', [off_set(ii)-1, p_centr(ceil(length(p_centr)/2)), 2, p_centr(end)- p_centr(ceil(length(p_centr)/2))], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii == 3 || ii ==4
        rectangle('position', [off_set(ii)-1, p_centr(1), 2, p_centr(ceil(length(p_centr)/2))- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    end
    
    hold on
    
    norm_val= max([nanmean(A_hist_pre(this_idx,:)) nanmean(A_hist_post(this_idx,:))]);
    plot((-nanmean(A_hist_pre(this_idx,:))./norm_val)+off_set(ii), p_centr,   'linewidth', 2, 'Color', f_ord(ii,:))
    plot((nanmean(A_hist_post(this_idx,:))./norm_val)+off_set(ii), p_centr,  'linewidth', 2, 'Color', [f_ord(ii,:) .5])
    
    if ii == 1
        %         text(off_set(ii)-.1, 97, '\leftarrow', 'HorizontalAlignment','right', 'Interpreter', 'TeX', 'fontsize', 8)
        text(off_set(ii)-.1, p_centr(end), 'pre', 'HorizontalAlignment','right', 'Interpreter', 'TeX', 'fontsize', 8, 'VerticalAlignment', 'baseline')
        %         text(off_set(ii)+.1, 97, '\rightarrow', 'HorizontalAlignment','left', 'Interpreter', 'TeX', 'fontsize', 8)
        text(off_set(ii)+.1, p_centr(end), 'post', 'HorizontalAlignment','left', 'Interpreter', 'TeX', 'fontsize', 8, 'VerticalAlignment', 'baseline')
        
    end
end
set(gca, 'xtick', off_set, 'ytick', [0 100])
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
ylabel('position (cm)')
axis('square')
title('REM assembly tuning')


subplot(1,2,2)
off_set = 0:2.5:10;
for ii = 1:5
    if ii == 1
        this_idx = lt1_idx;
    elseif ii == 2
        this_idx = lt5_idx;
    elseif ii == 3
        this_idx = H1_idx;
    elseif ii == 4
        this_idx = H5_idx;
    elseif ii == 5
        this_idx = HS_idx;
    end
    % add 'walls'
    if ii == 1 || ii == 2
        rectangle('position', [off_set(ii)-1, p_centr(1), 2, p_centr(end)- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii == 5
        rectangle('position', [off_set(ii)-1, p_centr(ceil(length(p_centr)/2)), 2, p_centr(end)- p_centr(ceil(length(p_centr)/2))], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    elseif ii == 3 || ii ==4
        rectangle('position', [off_set(ii)-1, p_centr(1), 2, p_centr(ceil(length(p_centr)/2))- p_centr(1)], 'FaceColor', [.9 .9 .9 .9], 'EdgeColor', [1 1 1]);
    end
    
    hold on
    
    norm_val= max([nanmean(n_pre_hist(this_idx,:)) nanmean(n_post_hist(this_idx,:))]);
    plot((-nanmean(n_pre_hist(this_idx,:))./norm_val)+off_set(ii), p_centr,   'linewidth', 2, 'Color', f_ord(ii,:))
    plot((nanmean(n_post_hist(this_idx,:))./norm_val)+off_set(ii), p_centr,  'linewidth', 2, 'Color', [f_ord(ii,:) .5])
    
    if ii == 1
        %         text(off_set(ii)-.1, 97, '\leftarrow', 'HorizontalAlignment','right', 'Interpreter', 'TeX', 'fontsize', 8)
        text(off_set(ii)-.1, p_centr(end), 'pre', 'HorizontalAlignment','right', 'Interpreter', 'TeX', 'fontsize', 8, 'VerticalAlignment', 'baseline')
        %         text(off_set(ii)+.1, 97, '\rightarrow', 'HorizontalAlignment','left', 'Interpreter', 'TeX', 'fontsize', 8)
        text(off_set(ii)+.1, p_centr(end), 'post', 'HorizontalAlignment','left', 'Interpreter', 'TeX', 'fontsize', 8, 'VerticalAlignment', 'baseline')
        
    end
end
set(gca, 'xtick', off_set, 'ytick', [0 100])
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
ylabel('position (cm)')
axis('square')
title({'REM assembly tuning'; 'wake normalized'})


%
% figure(99)
% clf
set(gcf, 'Units', 'centimeters', 'Position', [0 0  30 20])

cfg_set.ft_size = 12;
SetFigure(cfg_set, gcf)


%% Fig 5 example

%% plot a specific example
S_idx = []; 
for ii = length(A_out):-1:1
    if strcmpi(A_out{ii}{1}.info.subject, 'pv1069') &&  strcmpi(A_out{ii}{1}.info.session, 'LTD5')
        S_idx(ii) = true;
    else
        S_idx(ii) = false;
    end
    fprintf('%s   - %s\n', A_out{ii}{1}.info.subject, A_out{ii}{1}.info.session)
end

S_idx = find(S_idx);


% figure(9910)
% clf

    this_data = A_out{S_idx}{1};
    
    p_cent = []; p_rate = [];
    
    for ii = length(this_data.P_pos):-1:1
        p_cent(ii) = this_data.map{ii}.cent_z;
        p_rate(ii) = length(this_data.map{ii}.MI);%length(this_data.P_loc{ii}.loc_time);
    end
%     
    %p_rate(p_cent > -1.95) = 0;
    [~, s_idx] = sort(p_rate, 'descend');
    p_rank = s_idx;
    
    if length(this_data.P_pos) >3
        p_idx = p_rank(1:4);
    else
        p_idx = p_rank(1:length(this_data.P_pos));
    end
%     
%     c_ord = [67, 127, 151; 132 147 35; 255 179 13; 253 22 26]/255;
%     
%     
%     
%         ax1(1) = subplot(5,1,1);
%         cla; hold on
%                 rectangle('position', [0 0, this_data.behav.time(end)/1000, 50], 'Facecolor',  [c_ord(4,:) .6], 'Edgecolor',  [c_ord(4,:) 0])
%                 rectangle('position', [0 50, this_data.behav.time(end)/1000, 50], 'Facecolor',  [c_ord(1,:) .6], 'Edgecolor',  [c_ord(1,:) 0])
% 
%         plot(this_data.behav.time/1000, this_data.behav.position(:,1), 'k', 'linewidth', 3)
        
        
%                 ax1(2) = subplot(4,1,2:3);
r_idx = [1:4];%  [1 10];% [4 6 10]; %[1,3,4]; 
% r_idx = [1:length(p_rank)];
                c_ord_p =MS_linspecer(length(p_rank(r_idx)));

                % c_ord_p =MS_linspecer(40);
%                 c_ord_p = c_ord_p(5:floor(40/length(r_idx)):40,:); 
%                 c_ord_p = [c_ord_p(5,:);  c_ord_p(35,:); c_ord_p(40,:)];
% c_ord_p = linspecer(2);
               ax =  MS_Asmbly_plot_raster_figure(this_data, [], p_rank(r_idx), c_ord_p);
                linkprop(ax, 'XLim');
                xlim([0 100])
                MS_Asmbly_plot_raster_ReAct_figure(this_data, [], 'REM_Post_data', p_rank(r_idx), c_ord_p)
                xlim([0 100])

                figure(8988); clf; 
                if length(r_idx) <5
                    mm = 2;
                    nn = 2;
                else
                mm = ceil(length(r_idx)/2); 
                                nn = ceil(length(r_idx)/2); 
                end
                for ii = 1:length(r_idx)
                    subplot(mm,nn,ii); 
                    imagesc(this_data.map{p_rank(r_idx(ii))}.map); 
                    ylabel({num2str(p_rank(r_idx(ii))) ; [num2str(this_data.REM_Post_stats.rate(p_rank(r_idx(ii)))) '  |  ' num2str(this_data.REM_Post_stats.rate_p(p_rank(r_idx(ii))))  ]}); 
                end


% make the stem plots

figure(30311)
subplot(2,2,1)
    cla
    hold on
    ii = p_rank(r_idx(1)); 
    stem(this_data.P_temp(:,ii), 'color', [.8 .8 .8 .2])
    stem(this_data.P_pos{ii}, this_data.P_temp(this_data.P_pos{ii},ii), 'color', c_ord(1,:), 'MarkerFaceColor', c_ord(1,:), 'LineWidth', 1.5)
    view(90,90)
    ylim([-.2 .5])

    subplot(2,2,2)
    cla
    hold on
    ii = p_rank(r_idx(2)); 
    stem(this_data.P_temp(:,ii), 'color', [.8 .8 .8 .2])
    stem(this_data.P_pos{ii}, this_data.P_temp(this_data.P_pos{ii},ii), 'color', c_ord(2,:), 'MarkerFaceColor', c_ord(2,:), 'LineWidth', 1.5)
    view(90,90)
    ylim([-.2 .5])

    subplot(2,2,3)
    cla
    hold on
    ii = p_rank(r_idx(3)); 
    stem(this_data.P_temp(:,ii), 'color', [.8 .8 .8 .2])
    stem(this_data.P_pos{ii}, this_data.P_temp(this_data.P_pos{ii},ii), 'color', c_ord(3,:), 'MarkerFaceColor', c_ord(3,:), 'LineWidth', 1.5)
    view(90,90)
        ylim([-.2 .5])


    subplot(2,2,4)
    cla
    hold on
    ii = p_rank(r_idx(4)); 
    stem(this_data.P_temp(:,ii), 'color', [.8 .8 .8 .2])
    stem(this_data.P_pos{ii}, this_data.P_temp(this_data.P_pos{ii},ii), 'color', c_ord(4,:), 'MarkerFaceColor', c_ord(4,:), 'LineWidth', 1.5)
    view(90,90)
    ylim([-.2 .5])






%         yline(this_data.REM_Post_stats.R_thresh, '--', 'color', [.7 .7 .7], 'linewidth', 0.3)
%     yline(log10(this_data.REM_Post_stats.R_thresh), '--', 'color', [.7 .7 .7], 'linewidth', 0.3)
    
%     hold on
%     for ii = length(p_idx):-1:1
%         this_proj = this_data.P_proj(p_idx(ii),:);
%         plot(this_data.wake_tvec, log10(this_proj), 'color', c_ord(ii,:), 'linewidth', .5)
%     end
%     
%     ylim([0 inf])
%     y_val = get(gca, 'YTick');
%     %         y_val = [0 1 2 2.3979];
%     % yTick', [0 1 2 2.3979]
%     set(gca, 'YTickLabel', 10.^y_val, 'linewidth', 1);
%     
%         ylabel({'assembly strength'})
%     set(gca, 'xtick', []);
%     
        
        
%     
%         ax2(iD) = subplot(2,1,2);
% 
%     %         yline(this_data.REM_Post_stats.R_thresh, '--', 'color', [.7 .7 .7], 'linewidth', 0.3)
%     yline(log10(this_data.REM_Post_stats.R_thresh), '--', 'color', [.7 .7 .7], 'linewidth', 0.3)
%     
%     hold on
%     for ii = length(p_idx):-1:1
%         this_proj = this_data.REM_Post_proj(p_idx(ii),:);
%         plot(this_data.REM_Post_tvec, log10(this_proj), 'color', c_ord(ii,:), 'linewidth', .5)
%     end
%     
%     ylim([0 inf])
%     y_val = get(gca, 'YTick');
%     %         y_val = [0 1 2 2.3979];
%     % yTick', [0 1 2 2.3979]
%     set(gca, 'YTickLabel', 10.^y_val, 'linewidth', 1);
%     
%         ylabel({'assembly strength'})
%     set(gca, 'xtick', []);
%     
%         xlim([0 60]);
        
        
  
%     for ii = length(p_idx):-1:1
%         this_proj_idx = find(this_data.REM_Post_proj(p_idx(ii),:) > this_data.REM_Post_stats.R_thresh);
%         line([this_data.REM_Post_tvec(this_proj_idx); this_data.REM_Post_tvec(this_proj_idx)] , [ii-.5; ii+.5], 'color', c_ord(ii,:), 'linewidth', .5)
%     end
%     
%     set(gca, 'ytick', [])
%     
%     if iD == 1
%         xlim([0 60]);
%         xlabel('time (s)')
%         set(gca, 'xtick', [0 60])
   
    
% linkaxes(ax, 'y')
% set(gcf,'PaperUnits','inches', 'Units', 'inches');
% set(gcf, 'position', [5 5 7 2])
% set(gcf,'PaperSize', [7, 2]);
% set(gca,'xlimmode','manual','ylimmode','manual')
% print(gcf, '-dpdf', [fig_dir filesep 'Fig3_example.pdf'])

% save an h5 with the data for plotting.
% fname = ['Fig3_assembly_.h5'];

% if exist(fname, 'file')
%     delete(fname)
% end

%% create a matrix of number of events and rate

% f_list(ii).name,bin_size, move_thresh, method, 256)

Num_sig_A_mat = NaN(3,3); 

for ii = 1:length(A_out)

A_temp = []; A_prog = []; wake_data = []; wake_tvec = [];
for iB = length(A_out{ii}{1}.bins):-1:1
    
    [A_temp{iB}, A_proj{iB}, wake_data{iB}, wake_tvec{iB}, A_opts{iB}] = MS_PCA_ICA_only(ms_trk_cut, move_idx, bin_s(iB),method, opts);
end


for iB = length(A_temp):-1:1
    
    [P_temp{iB}, P_proj{iB}, P_pos{iB}] = MS_Asmbly_select(A_temp{iB}, A_proj{iB}, 2);
    
    fprintf('[%.0f/%.0f = %.0f%%] Assemblies had cells with positive weights (%.2fs binsize)\n',size(P_temp{iB},2),size(A_temp{iB},2),  (size(P_temp{iB},2)/size(A_temp{iB},2))*100, bin_s(iB))
end




end

Num_sig_A_mat(1,1,ii) = size(A_out{ii}{1}.REM_Pre_proj,1); 
Num_sig_A_mat(1,1,ii) = size(A_out{ii}{1}.REM_Pre_proj,1); 
Num_sig_A_mat(1,1,ii) = size(A_out{ii}{1}.REM_Pre_proj,1); 
Num_sig_A_mat(1,1,ii) = size(A_out{ii}{1}.REM_Pre_proj,1); 




%% count number of PRE assemblies and percent sig reactivatied
n_p_a = []; n_p_r =[]; k_i = []; 
l_pre_r =[]; l_post_r = []; 
for ii = size(B_out,2):-1:1
    if strcmp(B_out{ii}{1}.info.session, 'LTD1')
        k_i(ii) =true; 
    else 
        k_i(ii) = false; 
    end

    l_pre_r(ii) = B_out{ii}{1}.REM_Pre_tvec(end) - B_out{ii}{1}.REM_Pre_tvec(1); 
    l_post_r(ii) = B_out{ii}{1}.REM_Post_tvec(end) - B_out{ii}{1}.REM_Post_tvec(1); 

    if isempty(B_out{ii}{1}.pREM_Wake_proj)
        n_p_a(ii) =0;
        n_p_r(ii) =0;

    else
        n_p_a(ii) = size(B_out{ii}{1}.pREM_proj,1);
        n_p_r(ii) = sum(sum(B_out{ii}{1}.pREM_Wake_proj > B_out{ii}{1}.pREM_stats.R_thresh)>0);

    end
        
    s_t{ii} = B_out{ii}{1}.info.session; 

end

fprintf('Pre REM dur: %.2f +/- %.2f\n', mean(l_pre_r), std(l_pre_r))
fprintf('Post REM dur: %.2f +/- %.2f\n', mean(l_post_r), std(l_post_r))
fprintf('Pre REM prct: %.2f +/- %.2f\n', mean(l_pre_r/120), std(l_pre_r/120))
fprintf('Post REM prct: %.2f +/- %.2f\n', mean(l_post_r/120), std(l_post_r/120))

%% plot a specific example
S_idx = []; 
for ii = length(B_out):-1:1
    if strcmpi(B_out{ii}{1}.info.subject, 'pv1060') &&  strcmpi(B_out{ii}{1}.info.session, 'LTD1')
        S_idx(ii) = true;
    else
        S_idx(ii) = false;
    end
    fprintf('%s   - %s\n', B_out{ii}{1}.info.subject, B_out{ii}{1}.info.session)
end

S_idx = find(S_idx);


% figure(9910)
% clf

    this_data = B_out{S_idx}{1};
    
    p_cent = []; p_rate = [];
    
    for ii = size(this_data.pREM_temp,2):-1:1
        p_cent(ii) = this_data.map{ii}.cent_z;
        p_rate(ii) = length(this_data.P_loc{ii}.loc_time);
    end
%     
    p_rate(p_cent > -1.95) = 0;
    [~, s_idx] = sort(p_rate, 'descend');
    p_rank = s_idx;
    
    if length(this_data.P_pos) >3
        p_idx = p_rank(1:4);
    else
        p_idx = p_rank(1:length(this_data.P_pos));
    end

        
r_idx = 1:4;%  [1 10];% [4 6 10]; %[1,3,4]; 
% r_idx = [1:length(p_rank)];
                c_ord_p =MS_linspecer(length(p_rank(r_idx)));

                % c_ord_p =MS_linspecer(40);
%                 c_ord_p = c_ord_p(5:floor(40/length(r_idx)):40,:); 
%                 c_ord_p = [c_ord_p(5,:);  c_ord_p(35,:); c_ord_p(40,:)];
% c_ord_p = linspecer(2);
                MS_Asmbly_plot_raster_figure(this_data, [], p_rank(r_idx), c_ord_p)
                
                MS_Asmbly_plot_raster_ReAct_figure(this_data, [], 'REM_Post_data', p_rank(r_idx), c_ord_p)

                figure(8988); clf; 
                if length(r_idx) <5
                    mm = 2;
                    nn = 2;
                else
                mm = ceil(length(r_idx)/4); 
                                nn = ceil(length(r_idx)/4); 
                end
                for ii = 1:length(r_idx)
                    subplot(mm,nn,ii); 
                    imagesc(this_data.map{p_rank(r_idx(ii))}.map); 
                    ylabel({num2str(p_rank(r_idx(ii))) ; [num2str(this_data.REM_Post_stats.rate(p_rank(r_idx(ii)))) '  |  ' num2str(this_data.REM_Post_stats.rate_p(p_rank(r_idx(ii))))  ]}); 
                end