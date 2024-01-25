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
    Pipline_Asmbly_plot(A_out{ii}, fig_dir);
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
    Pipline_Asmbly_plot(J_out{ii}, fig_dir);
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
bar([mean(Pre_n_Asmbly(n_idx,iB)), mean(Post_n_Asmbly(n_idx,iB))], 'FaceColor', c_ord(1,:), 'EdgeColor', c_ord(1,:));
hold on
eb = errorbar([mean(Pre_n_Asmbly(n_idx,iB)), mean(Post_n_Asmbly(n_idx,iB))], [MS_SEM(Pre_n_Asmbly(n_idx,iB)) ,MS_SEM(Post_n_Asmbly(n_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title(['Novel (' num2str(bin_size(iB)) 's bins)'])
ylabel('# assemblies')
ylim([0 max_n_A(iB)])


subplot(n,m,2)
cla
bar([mean(Pre_n_Asmbly(~n_idx,iB)), mean(Post_n_Asmbly(~n_idx,iB))], 'FaceColor', c_ord(2,:), 'EdgeColor', c_ord(2,:));
hold on
eb = errorbar([mean(Pre_n_Asmbly(~n_idx,iB)), mean(Post_n_Asmbly(~n_idx,iB))], [MS_SEM(Pre_n_Asmbly(~n_idx,iB)) ,MS_SEM(Post_n_Asmbly(~n_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Fam')
ylim([0 max_n_A(iB)])


subplot(n,m,3)
cla
bar([mean(Pre_n_Asmbly(~a_idx,iB)), mean(Post_n_Asmbly(~a_idx,iB))], 'FaceColor', c_ord(3,:), 'EdgeColor', c_ord(3,:));
hold on
eb = errorbar([mean(Pre_n_Asmbly(~a_idx,iB)), mean(Post_n_Asmbly(~a_idx,iB))], [MS_SEM(Pre_n_Asmbly(~a_idx,iB)) ,MS_SEM(Post_n_Asmbly(~a_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Linear')
ylim([0 max_n_A(iB)])

subplot(n,m,4)
cla
bar([mean(Pre_n_Asmbly(a_idx,iB)), mean(Post_n_Asmbly(a_idx,iB))], 'FaceColor', c_ord(4,:), 'EdgeColor', c_ord(4,:));
hold on
eb = errorbar([mean(Pre_n_Asmbly(a_idx,iB)), mean(Post_n_Asmbly(a_idx,iB))], [MS_SEM(Pre_n_Asmbly(a_idx,iB)) ,MS_SEM(Post_n_Asmbly(a_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Anxiety')
ylim([0 max_n_A(iB)])

subplot(n,m,5)
cla
bar([mean(Pre_n_Asmbly(HS_idx,iB)), mean(Post_n_Asmbly(HS_idx,iB))], 'FaceColor', c_ord(5,:), 'EdgeColor', c_ord(5,:));
hold on
eb = errorbar([mean(Pre_n_Asmbly(HS_idx,iB)), mean(Post_n_Asmbly(HS_idx,iB))], [MS_SEM(Pre_n_Asmbly(HS_idx,iB)) ,MS_SEM(Post_n_Asmbly(HS_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('A. Switch')
ylim([0 max_n_A(iB)])


% same for the rate

subplot(n,m,m+1)
cla
bar([mean(Pre_r_Asmbly(n_idx,iB)), mean(Post_r_Asmbly(n_idx,iB))], 'FaceColor', c_ord(1,:), 'EdgeColor', c_ord(1,:));
hold on
eb = errorbar([mean(Pre_r_Asmbly(n_idx,iB)), mean(Post_r_Asmbly(n_idx,iB))], [MS_SEM(Pre_r_Asmbly(n_idx,iB)) ,MS_SEM(Post_r_Asmbly(n_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title(['Novel (' num2str(bin_size(iB)) 's bins)'])
ylabel('ReAct/min')
ylim([0 max_r_A(iB)])



subplot(n,m,m+2)
cla
bar([mean(Pre_r_Asmbly(~n_idx,iB)), mean(Post_r_Asmbly(~n_idx,iB))], 'FaceColor', c_ord(2,:), 'EdgeColor', c_ord(2,:));
hold on
eb = errorbar([mean(Pre_r_Asmbly(~n_idx,iB)), mean(Post_r_Asmbly(~n_idx,iB))], [MS_SEM(Pre_r_Asmbly(~n_idx,iB)) ,MS_SEM(Post_r_Asmbly(~n_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Fam')
ylim([0 max_r_A(iB)])


subplot(n,m,m+3)
cla
bar([mean(Pre_r_Asmbly(~a_idx,iB)), mean(Post_r_Asmbly(~a_idx,iB))], 'FaceColor', c_ord(3,:), 'EdgeColor', c_ord(3,:));
hold on
eb = errorbar([mean(Pre_r_Asmbly(~a_idx,iB)), mean(Post_r_Asmbly(~a_idx,iB))], [MS_SEM(Pre_r_Asmbly(~a_idx,iB)) ,MS_SEM(Post_r_Asmbly(~a_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Linear')
ylim([0 max_r_A(iB)])


subplot(n,m,m+4)
cla
bar([mean(Pre_r_Asmbly(a_idx,iB)), mean(Post_r_Asmbly(a_idx,iB))], 'FaceColor', c_ord(4,:), 'EdgeColor', c_ord(4,:));
hold on
eb = errorbar([mean(Pre_r_Asmbly(a_idx,iB)), mean(Post_r_Asmbly(a_idx,iB))], [MS_SEM(Pre_r_Asmbly(a_idx,iB)) ,MS_SEM(Post_r_Asmbly(a_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Anxiety')
ylim([0 max_r_A(iB)])


subplot(n,m,m+5)
cla
bar([mean(Pre_r_Asmbly(HS_idx,iB)), mean(Post_r_Asmbly(HS_idx,iB))], 'FaceColor', c_ord(5,:), 'EdgeColor', c_ord(5,:));
hold on
eb = errorbar([mean(Pre_r_Asmbly(HS_idx,iB)), mean(Post_r_Asmbly(HS_idx,iB))], [MS_SEM(Pre_r_Asmbly(HS_idx,iB)) ,MS_SEM(Post_r_Asmbly(HS_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('A. Switch')
ylim([0 max_r_A(iB)])



% repeat for assemblies with coherent spatial tuning. 


subplot(n,m,(m*2)+1)
cla
bar([nanmean(S_Pre_n_Asmbly(n_idx,iB)), nanmean(S_Post_n_Asmbly(n_idx,iB))], 'FaceColor', c_ord_s(1,:), 'EdgeColor', c_ord_s(1,:));
hold on
eb = errorbar([nanmean(S_Pre_n_Asmbly(n_idx,iB)), nanmean(S_Post_n_Asmbly(n_idx,iB))], [MS_SEM(S_Pre_n_Asmbly(n_idx,iB)) ,MS_SEM(S_Post_n_Asmbly(n_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title(['Novel (' num2str(bin_size(iB)) 's bins)'])
ylabel('# assemblies')
ylim([0 max_n_SA(iB)])



subplot(n,m,(m*2)+2)
cla
bar([nanmean(S_Pre_n_Asmbly(~n_idx,iB)), nanmean(S_Post_n_Asmbly(~n_idx,iB))], 'FaceColor', c_ord_s(2,:), 'EdgeColor', c_ord_s(2,:));
hold on
eb = errorbar([nanmean(S_Pre_n_Asmbly(~n_idx,iB)), nanmean(S_Post_n_Asmbly(~n_idx,iB))], [MS_SEM(S_Pre_n_Asmbly(~n_idx,iB)) ,MS_SEM(S_Post_n_Asmbly(~n_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Fam')
ylim([0 max_n_SA(iB)])


subplot(n,m,(m*2)+3)
cla
bar([nanmean(S_Pre_n_Asmbly(~a_idx,iB)), nanmean(S_Post_n_Asmbly(~a_idx,iB))], 'FaceColor', c_ord_s(3,:), 'EdgeColor', c_ord_s(3,:));
hold on
eb = errorbar([nanmean(S_Pre_n_Asmbly(~a_idx,iB)), nanmean(S_Post_n_Asmbly(~a_idx,iB))], [MS_SEM(S_Pre_n_Asmbly(~a_idx,iB)) ,MS_SEM(S_Post_n_Asmbly(~a_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Linear')
ylim([0 max_n_SA(iB)])

subplot(n,m,(m*2)+4)
cla
bar([nanmean(S_Pre_n_Asmbly(a_idx,iB)), nanmean(S_Post_n_Asmbly(a_idx,iB))], 'FaceColor', c_ord_s(4,:), 'EdgeColor', c_ord_s(4,:));
hold on
eb = errorbar([nanmean(S_Pre_n_Asmbly(a_idx,iB)), nanmean(S_Post_n_Asmbly(a_idx,iB))], [MS_SEM(S_Pre_n_Asmbly(a_idx,iB)) ,MS_SEM(S_Post_n_Asmbly(a_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Anxiety')
ylim([0 max_n_SA(iB)])


subplot(n,m,(m*2)+5)
cla
bar([nanmean(S_Pre_n_Asmbly(HS_idx,iB)), nanmean(S_Post_n_Asmbly(HS_idx,iB))], 'FaceColor', c_ord_s(5,:), 'EdgeColor', c_ord_s(5,:));
hold on
eb = errorbar([nanmean(S_Pre_n_Asmbly(HS_idx,iB)), nanmean(S_Post_n_Asmbly(HS_idx,iB))], [MS_SEM(S_Pre_n_Asmbly(HS_idx,iB)) ,MS_SEM(S_Post_n_Asmbly(HS_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('A. Switch')
ylim([0 max_n_SA(iB)])


% same for the rate

subplot(n,m,(m*3)+1)
cla
bar([nanmean(S_Pre_r_Asmbly(n_idx,iB)), nanmean(S_Post_r_Asmbly(n_idx,iB))], 'FaceColor', c_ord_s(1,:), 'EdgeColor', c_ord_s(1,:));
hold on
eb = errorbar([nanmean(S_Pre_r_Asmbly(n_idx,iB)), nanmean(S_Post_r_Asmbly(n_idx,iB))], [MS_SEM(S_Pre_r_Asmbly(n_idx,iB)) ,MS_SEM(S_Post_r_Asmbly(n_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title(['Novel (' num2str(bin_size(iB)) 's bins)'])
ylabel('ReAct/min')
ylim([0 max_r_SA(iB)])



subplot(n,m,(m*3)+2)
cla
bar([nanmean(S_Pre_r_Asmbly(~n_idx,iB)), nanmean(S_Post_r_Asmbly(~n_idx,iB))], 'FaceColor', c_ord_s(2,:), 'EdgeColor', c_ord_s(2,:));
hold on
eb = errorbar([nanmean(S_Pre_r_Asmbly(~n_idx,iB)), nanmean(S_Post_r_Asmbly(~n_idx,iB))], [MS_SEM(S_Pre_r_Asmbly(~n_idx,iB)) ,MS_SEM(S_Post_r_Asmbly(~n_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Fam')
ylim([0 max_r_SA(iB)])


subplot(n,m,(m*3)+3)
cla
bar([nanmean(S_Pre_r_Asmbly(~a_idx,iB)), nanmean(S_Post_r_Asmbly(~a_idx,iB))], 'FaceColor', c_ord_s(3,:), 'EdgeColor', c_ord_s(3,:));
hold on
eb = errorbar([nanmean(S_Pre_r_Asmbly(~a_idx,iB)), nanmean(S_Post_r_Asmbly(~a_idx,iB))], [MS_SEM(S_Pre_r_Asmbly(~a_idx,iB)) ,MS_SEM(S_Post_r_Asmbly(~a_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Linear')
ylim([0 max_r_SA(iB)])


subplot(n,m,(m*3)+4)
cla
bar([nanmean(S_Pre_r_Asmbly(a_idx,iB)), nanmean(S_Post_r_Asmbly(a_idx,iB))], 'FaceColor', c_ord_s(4,:), 'EdgeColor', c_ord_s(4,:));
hold on
eb = errorbar([nanmean(S_Pre_r_Asmbly(a_idx,iB)), nanmean(S_Post_r_Asmbly(a_idx,iB))], [MS_SEM(S_Pre_r_Asmbly(a_idx,iB)) ,MS_SEM(S_Post_r_Asmbly(a_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Anxiety')
ylim([0 max_r_SA(iB)])


subplot(n,m,(m*3)+5)
cla
bar([nanmean(S_Pre_r_Asmbly(HS_idx,iB)), nanmean(S_Post_r_Asmbly(HS_idx,iB))], 'FaceColor', c_ord_s(5,:), 'EdgeColor', c_ord_s(5,:));
hold on
eb = errorbar([nanmean(S_Pre_r_Asmbly(HS_idx,iB)), nanmean(S_Post_r_Asmbly(HS_idx,iB))], [MS_SEM(S_Pre_r_Asmbly(HS_idx,iB)) ,MS_SEM(S_Post_r_Asmbly(HS_idx,iB))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('A. Switch')
ylim([0 max_r_SA(iB)])


end