%% MASTER_Asmbly

if strcmp(computer, 'GLNXA64')
    
    codebase_dir = '/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared';
    ca_dir = '/home/williamslab/Documents/Github/CEH2';
    oasis_dir = '/home/williamslab/Documents/Github/OASIS_matlab';
    
    code_dir = '/home/williamslab/Documents/Github/Dos-Santos Assembly ICA/Dos-Santos Assembly ICA';
    
    RnR_dir = '/home/williamslab/Documents/Github/RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    main_dir = '/home/williamslab';
    
    
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
bin_size = [.5];
%%
%
% cd('C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\PC9')
% fig_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\PC9\checks';
%
% % cd('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter')
%
% f_list = dir('*data*');
%

%
% A_out = [];
% session = []; novel_idx = []; anx_idx = []; HS_idx = [];
% method = 'grosmark';
%
%
%
%
% for ii = 13:length(f_list)
%     session{ii} = f_list(ii).name;
%
%     % compute assemblies and related ReActs
%     A_out{ii} = Pipeline_Asmbly(f_list(ii).name,bin_size, move_thresh, method);
%
%     % Summary plots
%             Pipline_Asmbly_plot(A_out{ii}, [fig_dir filesep method]);
%     close all
%
%     if ~isempty(strfind(f_list(ii).name, 'HATDS'))
%         HS_idx(ii) = 1;
%     else
%         HS_idx(ii) = 0;
%     end
%
%     if ~isempty(strfind(f_list(ii).name, 'D1')) %|| ~isempty(strfind(f_list(ii).name, 'HATDS'))
%         novel_idx(ii) = 1;
%     end
%
%     if ~isempty(strfind(f_list(ii).name, 'D5'))
%         novel_idx(ii) = 0;
%     end
%
%     if ~isempty(strfind(f_list(ii).name, 'HAT'))
%         anx_idx(ii) = 1;
%     else
%         anx_idx(ii) = 0;
%     end
%
% end
%
% A_all = A_out;
%
% save(['C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\inter\A_out_' method '.mat'], 'A_out')
%
% %%  Run again but for EVV data
% data_dir = ('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV');
% fig_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV\checks';
% % cd('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter')
% cd(data_dir);
% j20_list = dir('*data*');
%
% method = 'grosmark';
% fig_dir = ['C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV\checks\' method];
%
% J_out = [];
% J20_session = []; J20_novel_idx = []; D3_idx = [];
% for ii = 1:length(j20_list)
%     J20_session{ii} = j20_list(ii).name;
%     cd(data_dir)
%
%     J_out{ii} = Pipeline_Asmbly(j20_list(ii).name,bin_size, move_thresh, method);
%
%     % Summary plots
%     %         Pipline_Asmbly_plot(J_out{ii}, fig_dir);
%
%
%     close all
%
%     if ~isempty(strfind(j20_list(ii).name, 'D1')) %|| ~isempty(strfind(f_list(ii).name, 'HATDS'))
%         J20_novel_idx(ii) = 1;
%         D3_idx(ii) = 0;
%     end
%
%     if ~isempty(strfind(j20_list(ii).name, 'D3'))
%         D3_idx(ii) = 1;
%         J20_novel_idx(ii) = 0;
%     end
%
%     if ~isempty(strfind(j20_list(ii).name, 'D5'))
%         J20_novel_idx(ii) = 0;
%         D3_idx(ii) = 0;
%     end
%
% end
%
% J_all = J_out;
% % A_out = J_out;
% save(['C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\inter\J_out_' method '.mat'], 'J_out')
%%  Run again but for EVV data
data_dir = ('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV');
fig_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV\checks';
% cd('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter')
cd(data_dir);
j20_list = dir('*data.mat*');

for ii = length(j20_list):-1:1
   if contains(j20_list(ii).name, 'LTD3')
       keep_idx(ii) = 0; 
   else
       keep_idx(ii) = 1; 
   end
    
end

j20_list(~keep_idx) = [];  

method = 'binary';

J_out = []; JP_out = [];
J20_session = []; J20_novel_idx = []; D3_idx = []; D5_idx = [];
for ii = 1:length(j20_list)
    
    
    J20_session{ii} = j20_list(ii).name;
    cd(data_dir)
        
        %             J_out{ii} = Pipeline_Asmbly(j20_list(ii).name,bin_size, move_thresh, method);
        %         JP_out{ii} = Pipeline_Asmbly_place(j20_list(ii).name,bin_size, move_thresh, method);
        
%         J_out{ii} = Pipeline_Asmbly_top_cells(j20_list(ii).name,bin_size, move_thresh, method);
%         J_out{ii} = Pipeline_Asmbly_append_preA(J_out{ii});
        
%         Pipline_Asmbly_plot(J_out{ii}, [fig_dir filesep method filesep 'best']);
        
        
        % Summary plots
%         Pipline_Asmbly_plot(J_out{ii}, [fig_dir filesep method]);
        close all

        % Pipline_Asmbly_plot(JP_out{ii}, [fig_dir filesep method filesep 'place']);
    % close all
    
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
J20_novel_idx = logical(J20_novel_idx);
D3_idx = logical(D3_idx);
D5_idx = logical(~D3_idx & ~J20_novel_idx);


if ~isempty(J_out)
    save(['C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\inter\J_out_' method '.mat'], 'J_out')
end
if ~isempty(JP_out)
    save(['C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\inter\JP_out_' method '.mat'], 'JP_out')
end
%%
inter_dir = strrep([main_dir 'Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\PC9'], '\', filesep);
fig_dir = [inter_dir filesep 'checks'];
% cd('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter')

cd(inter_dir);

f_list = dir('*data*');

A_out = []; P_out = [];
session = []; novel_idx = []; anx_idx = []; HS_idx = [];
method = 'binary';

for ii = 1:length(f_list)
    session{ii} = f_list(ii).name;
    % compute assemblies and related ReActs
    %         A_out{ii} = Pipeline_Asmbly(f_list(ii).name,bin_size, move_thresh, method);
    %     P_out{ii} = Pipeline_Asmbly_place(f_list(ii).name,bin_size, move_thresh, method);
    %
    B_out{ii} = Pipeline_Asmbly_top_cells(f_list(ii).name,bin_size, move_thresh, method);
    %
    %         B_out{ii} = Pipeline_Asmbly_append_SWS(f_list(ii).name, B_out{ii});
    %
    B_out{ii} = Pipeline_Asmbly_append_preA(B_out{ii});
    %
    % % Summary plots
    % %                 Pipline_Asmbly_plot(A_out{ii}, [fig_dir filesep method]);
    % %     Pipline_Asmbly_plot(P_out{ii}, [fig_dir filesep method filesep 'place']);
    %
    Pipline_Asmbly_plot(B_out{ii}, [fig_dir filesep method filesep 'best']);
    %
    % Pipline_Asmbly_plot_SWS(B_out{ii}, [fig_dir filesep method filesep 'best_SWS']);
    
    
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
        B_out{ii} = Pipeline_Asmbly_append_preA(B_out{ii});
        
    end
    
end
novel_idx = logical(novel_idx);
anx_idx = logical(anx_idx);
HS_idx = logical(HS_idx);

lt1_idx = novel_idx & ~anx_idx & ~HS_idx;
lt5_idx = ~novel_idx & ~anx_idx & ~HS_idx;

H1_idx = novel_idx & anx_idx & ~HS_idx;
H5_idx = ~novel_idx & anx_idx & ~HS_idx;

% A_all = A_out;
% A_out = A_out(1:11);
% % if ~isempty(A_out)
% %     save(['C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\inter\A_out_' method '.mat'], 'A_out')
% % end
% % if ~isempty(P_out)
% save([main_dir  strrep('Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\inter\B_out_', '\', filesep) method '.mat'], 'B_out')
% % end
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

max_n_A = max([nanmean(Pre_n_Asmbly) ; nanmean(Post_n_Asmbly) ]);
max_n_A = max_n_A*1.5;
max_r_A = max([nanmean(Pre_r_Asmbly) ; nanmean(Post_r_Asmbly) ]);
max_r_A = max_r_A*1.5;

max_n_SA = max([nanmean(S_Pre_n_Asmbly) ; nanmean(S_Post_n_Asmbly) ]);
max_n_SA = max_n_SA*1.5;
max_r_SA = max([nanmean(S_Pre_r_Asmbly) ; nanmean(S_Post_r_Asmbly) ]);
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

J_n_idx = logical(J20_novel_idx);



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

fprintf('<strong> %s</strong> <strong> %0.02f</strong>\x00B1%0.2f wake assemblies | <strong> %0.02f</strong>\x00B1%0.2f pre | <strong> %0.02f</strong>\x00B1%0.2f post\n',...
    'overall', nanmean(wake_n_Asmbly(:, iB)), std(wake_n_Asmbly(:,iB)), nanmean(Pre_n_Asmbly(:, iB)), std(Pre_n_Asmbly(:,iB)), nanmean(Post_n_Asmbly(:, iB)), std(Post_n_Asmbly(:,iB)))

fprintf('<strong> %s</strong> <strong> %0.02f</strong>\x00B1%0.2f wake assemblies | <strong> %0.02f</strong>\x00B1%0.2f pre | <strong> %0.02f</strong>\x00B1%0.2f post\n',...
    'Novel', nanmean(wake_n_Asmbly(novel_idx, iB)), std(wake_n_Asmbly(novel_idx,iB)), nanmean(Pre_n_Asmbly(novel_idx, iB)), std(Pre_n_Asmbly(novel_idx,iB)), nanmean(Post_n_Asmbly(novel_idx, iB)), std(Post_n_Asmbly(novel_idx,iB)))

fprintf('<strong> %s</strong> <strong> %0.02f</strong>\x00B1%0.2f wake assemblies | <strong> %0.02f</strong>\x00B1%0.2f pre | <strong> %0.02f</strong>\x00B1%0.2f post\n',...
    'Familiar', nanmean(wake_n_Asmbly(~novel_idx, iB)), std(wake_n_Asmbly(~novel_idx,iB)), nanmean(Pre_n_Asmbly(~novel_idx, iB)), std(Pre_n_Asmbly(~novel_idx,iB)), nanmean(Post_n_Asmbly(~novel_idx, iB)), std(Post_n_Asmbly(~novel_idx,iB)))


fprintf('<strong> %s</strong> <strong> %0.02f</strong>\x00B1%0.2f wake assemblies | <strong> %0.02f</strong>\x00B1%0.2f pre | <strong> %0.02f</strong>\x00B1%0.2f post\n',...
    'LT', nanmean(wake_n_Asmbly(~anx_idx, iB)), std(wake_n_Asmbly(~anx_idx,iB)), nanmean(Pre_n_Asmbly(~anx_idx, iB)), std(Pre_n_Asmbly(~anx_idx,iB)), nanmean(Post_n_Asmbly(~anx_idx, iB)), std(Post_n_Asmbly(~anx_idx,iB)))

fprintf('<strong> %s</strong> <strong> %0.02f</strong>\x00B1%0.2f wake assemblies | <strong> %0.02f</strong>\x00B1%0.2f pre | <strong> %0.02f</strong>\x00B1%0.2f post\n',...
    'HAT', nanmean(wake_n_Asmbly(anx_idx, iB)), std(wake_n_Asmbly(anx_idx,iB)), nanmean(Pre_n_Asmbly(anx_idx, iB)), std(Pre_n_Asmbly(anx_idx,iB)), nanmean(Post_n_Asmbly(anx_idx, iB)), std(Post_n_Asmbly(anx_idx,iB)))


%% collect j20 and control reactivation strength

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


for iB = length(bin_size):-1:1
    for iA = size(J_out,2):-1:1
        
        %         wake
        J_wake_n(iA,iB) = size(J_out{iA}{iB}.P_proj,1);
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


fprintf('<strong> %s</strong> <strong> %0.02f</strong>\x00B1%0.2f wake assemblies | <strong> %0.02f</strong>\x00B1%0.2f pre | <strong> %0.02f</strong>\x00B1%0.2f post\n',...
    'J20 overall', nanmean(J_wake_n(:, iB)), std(J_wake_n(:,iB)), nanmean(J_Pre_n(:, iB)), std(J_Pre_n(:,iB)), nanmean(J_Post_a(:, iB)), std(J_Post_a(:,iB)))

fprintf('<strong> %s</strong> <strong> %0.02f</strong>\x00B1%0.2f wake assemblies | <strong> %0.02f</strong>\x00B1%0.2f pre | <strong> %0.02f</strong>\x00B1%0.2f post\n',...
    'J20 Novel', nanmean(J_wake_n(J20_novel_idx, iB)), std(J_wake_n(J20_novel_idx,iB)), nanmean(J_Pre_n(J20_novel_idx, iB)), std(J_Pre_n(J20_novel_idx,iB)), nanmean(J_Post_a(J20_novel_idx, iB)), std(J_Post_a(J20_novel_idx,iB)))

fprintf('<strong> %s</strong> <strong> %0.02f</strong>\x00B1%0.2f wake assemblies | <strong> %0.02f</strong>\x00B1%0.2f pre | <strong> %0.02f</strong>\x00B1%0.2f post\n',...
    'J20 D5', nanmean(J_wake_n(D5_idx, iB)), std(J_wake_n(D5_idx,iB)), nanmean(J_Pre_n(D5_idx, iB)), std(J_Pre_n(D5_idx,iB)), nanmean(J_Post_a(D5_idx, iB)), std(J_Post_a(D5_idx,iB)))

%% plot the J20 vs control novel VS familiar
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


J_n_idx = logical(J20_novel_idx);
J_f_idx = D3_idx | D5_idx;

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


% same but J20
subplot(2,4,4)
data_a = J_ReAct(J_n_idx, iB);
data_b = J_ReAct(D3_idx, iB);
data_c = J_ReAct(D5_idx, iB);
J_cord = summer(3);

hb = bar([nanmean(data_a), nanmean(data_b),nanmean(data_c)]', 'FaceColor', J_cord(1,:), 'EdgeColor', 'k');
hb.FaceColor = 'flat';
hb.CData(2,:) = J_cord(2,:);
hb.CData(3,:) = J_cord(3,:);
hold on
eb = errorbar([nanmean(data_a), nanmean(data_b), nanmean(data_c)], [MS_SEM(data_a) ,MS_SEM(data_b), MS_SEM(data_c)]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:3, 'XTickLabel', {'LT1', 'LT3', 'LT5'}, 'XTickLabelRotation', 45)
ylim(y_lim);
title('J20')


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

%% Collect the change in significant CCF
pre_sig_cff = []; post_sig_cff = [];
J_pre_sig_cff = []; J_post_sig_cff = [];
for iB = length(bin_size):-1:1
    
    for iS = size(A_out,2):-1:1
        
        pre_sig_cff(iS, iB) = A_out{iS}{iB}.REM_Pre_psig_cff;
        post_sig_cff(iS, iB) = A_out{iS}{iB}.REM_Post_psig_cff;
    end
    
    for iJ = size(J_out,2):-1:1
        
        J_pre_sig_cff(iJ, iB) = J_out{iJ}{iB}.REM_Pre_psig_cff;
        J_post_sig_cff(iJ, iB) = J_out{iJ}{iB}.REM_Post_psig_cff;
    end
    
    figure(10+iB)
    %     c_ord = MS_linspecer(4);
    ax(1) = subplot(2,5,1);
    cla
    MS_bar_w_err(pre_sig_cff(n_idx,iB), post_sig_cff(n_idx,iB), c_ord(1,:));
    set(gca, 'XTickLabel', {'pre', 'post'})
    title('novel');
    %     title(['CCF sig (bin: ' num2str(bin_size(iB)) ')'])
    ylabel({'Cross correlation'; '% significant of wake'})
    
    ax(2) = subplot(2,5,2);
    cla
    MS_bar_w_err(pre_sig_cff(f_idx,iB), post_sig_cff(f_idx,iB), c_ord(2,:))
    set(gca, 'XTickLabel', {'pre', 'post'})
    title('familiar');
    
    
    ax(3) = subplot(2,5,3);
    cla
    MS_bar_w_err(pre_sig_cff(lt1_idx | lt5_idx,iB), post_sig_cff(lt1_idx | lt5_idx,iB), c_ord(3,:))
    set(gca, 'XTickLabel', {'pre', 'post'})
    title('Linear');
    
    
    ax(4) = subplot(2,5,4);
    cla
    MS_bar_w_err(pre_sig_cff(H1_idx | H5_idx,iB), post_sig_cff(H1_idx | H5_idx,iB), c_ord(4,:))
    set(gca, 'XTickLabel', {'pre', 'post'})
    title('Anxiety');
    
    ax(5) =  subplot(2,5,5);
    cla
    MS_bar_w_err(pre_sig_cff(H1_idx | H5_idx,iB), post_sig_cff(H1_idx | H5_idx,iB), c_ord(5,:))
    set(gca, 'XTickLabel', {'pre', 'post'})
    title('switch');
    
    linkprop(ax, 'ylim')
    ylim([0 40])
    %     Square_subplots;
    SetFigure([], gcf)
    %     title('J20')
    saveas(gcf,[fig_dir filesep 'cff_summary_p' num2str(strrep(bin_size(iB), '.', 'p')) '.png'])
    
end


%% qualitfy reactivation saptial biases

data_in = B_out;
Pre_hist = []; Post_hist = [];
A_hist_pre = []; A_hist_post = [];
J20_hist_pre = []; J20_hist_post = [];
for kk = 1
    if kk == 1
        data_in = B_out;
    elseif kk == 2
        data_in = J_out;
    end
    
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
            elseif kk == 2
                J20_hist_pre(ii,:) = Pre_hist{ii,iB};
                J20_hist_post(ii,:) = Post_hist{ii,iB};
            end
            
        end
        
    end
    
    if kk == 1
        A_cnt_prct = all_cnt_prct;
        A_Pre_cnt_prct = Pre_cnt_prct;
        A_Post_cnt_prct = Post_cnt_prct;
        
    elseif kk == 2
        J_cnt_prct = all_cnt_prct;
        J_Pre_cnt_prct = Pre_cnt_prct;
        J_Post_cnt_prct = Post_cnt_prct;
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
stem(p_centr, nanmean(A_hist_pre(HS_idx,:))./max(nanmean(A_hist_pre(HS_idx,:))), 'color', a_ord(2,:));
stem(p_centr, nanmean(A_hist_post(HS_idx,:))./max(nanmean(A_hist_post(HS_idx,:))), 'color', a_ord(1,:));
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
    Pre_REM_nA = size(A_out{ii}{1}.REM_temp,2);
    
    if Pre_REM_nA == 0
        Pre_REM_wAct = 0;
    else
        Pre_REM_wAct = sum(A_out{ii}{1}.REM_Wake_proj > 8, 2);
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

conds  = {'LTD1 ', 'LTD5 ', 'HATD1 ', 'HATD5 ', 'HATDSwitch '};
mice = {'pv1043', 'pv1060', 'pv1069', 'pv1191', 'pv1192', 'pv1252', 'pv1254'};
Cond = []; Sub = []; nWake_A = [];nWake_R =[];  nPre_A = []; nPre_Act = []; nPre_SWS_A = []; nPost_SWS_A = [];
PreA_str = []; PostA_str= [];
nPre_sig = []; nPost_sig = [];
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

subplot(2,4,5)

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
    
    hold on
    norm_val= max([nanmean(A_hist_pre(this_idx,:)) nanmean(A_hist_post(this_idx,:))]);
    plot((-nanmean(A_hist_pre(this_idx,:))./norm_val)+off_set(ii), p_centr,   'linewidth', 2, 'Color', f_ord(ii,:))
    plot((nanmean(A_hist_post(this_idx,:))./norm_val)+off_set(ii), p_centr,  'linewidth', 2, 'Color', [f_ord(ii,:) .5])
end
set(gca, 'xtick', off_set, 'ytick', [0 100])
set(gca, 'XTickLabel', {'Novel', 'Familiar', 'HAT_{nov}', 'HAT_{fam}', 'HAT_{sw}'}, 'XTickLabelRotation', 45)
ylabel('position (cm)')
axis('square')




subplot(2,4,1)

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

print(gcf,  [fig_dir filesep   'Stats_summery_' strrep(num2str(B_out{1}{1}.info.bin), '.', 'p') 's_bin.pdf'], '-dpdf','-r0')



%%
% try to get the pre and post REM assemblies in one plot
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

lme_wake_A = fitlme(tbl,'nWake_A~Cond+(1|Sub)');
fprintf('nWake_A: coef test p = %2.4f\n', coefTest(lme_wake_A))

lme_pre_sig = fitlme(tbl,'nPre_sig~Cond+(1|Sub)');
fprintf('nPre_Sig: coef test p = %2.4f\n', coefTest(lme_pre_sig))

lme_post_sig = fitlme(tbl,'nPost_sig~Cond+(1|Sub)');
fprintf('nPost_Sig: coef test p = %2.4f\n', coefTest(lme_post_sig))

lme_pre_A = fitlme(tbl,'nPre_A~Cond+(1|Sub)');
fprintf('nPre_A: coef test p = %2.4f\n', coefTest(lme_pre_A))

lme_ReAct = fitlme(tbl,'ReAct_str ~Cond+(1|Sub)');
fprintf('ReAct_Str: coef test p = %2.4f\n', coefTest(lme_ReAct))

lme_wake_R = fitlme(tbl,'nWake_R ~Cond+(1|Sub)');
fprintf('nWake_R: coef test p = %2.4f\n', coefTest(lme_wake_R))
% stats = anova1(lme_wake_A);
% multcompare(stats)

%% quick plot of the reactivtion strength by condition




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


%% sample plots for PCA ICA

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