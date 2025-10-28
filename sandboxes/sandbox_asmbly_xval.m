%% MS_Asmbly_Shuff_tests


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

% clear the paths vars after adding.
clearvars *_dir -except main_dir

cd(c_d)


move_thresh  = 9;
bin_size = .5;

inter_dir = strrep([main_dir 'Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\PC9'], '\', filesep);
fig_dir = [inter_dir filesep 'checks'];

%%

method = 'binary';

warning off
load([main_dir  strrep('Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\inter\B_out_', '\', filesep) method '.mat'], 'B_out')
warning on

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


%% run quick stats on the number of assemblies and the shuffle metrics
rng(123, 'twister'); % for reproducibility.

opts = [];
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 500;
opts.threshold.permutations_percentile= 95;
opts.threshold.number_of_permutations= 500;


% init vars for number of assemblies per session and the shuff stats
a_num = NaN(size(A_out));
shuff_num = a_num;
wake_nA = []; wake_pre_nA = []; wake_post_nA = [];
J_sig = [];
for iA = length(A_out):-1:1

    % % number of awake assemblies and Shuff
    % wake_nA(iA) = size(A_out{iA}{1}.P_proj,1);
    % wake_shuf_nA(iA) = size(A_out{iA}{1}.P_proj,1);
    %
    % % number of awake in Pre REM
    % wake_pre_nA(iA) = sum(A_out{iA}{1}.REM_Pre_stats.p_val < 0.05);
    %
    % % number of awake in Post REM
    % wake_post_nA(iA) = sum(A_out{iA}{1}.REM_Post_stats.p_val < 0.05);


    % current assembly wieghts;
    A_temp = A_out{iA}{1}.P_temp;

    % apply the weights to other subjects



    for jj = length(A_out):-1:1



        % if jj == iA
        J_data = A_out{jj}{1}.REM_post_in;
        A_data = A_out{jj}{1}.REM_post_in;

        % test if A_temps can be found in J-data
        % _proj = assembly_activity(A_temp,J_data');

        A_alt_proj = [];
        for tt = size(A_alt_proj, 1):-1:1
            A_alt_proj(tt) = sum(A_alt_proj(tt,:)>9);
        end

        J_sig(iA, jj) = sum(J_proj > 0)/length(J_proj); % get the percentage of assemblies exceeding the activation cut off.

        J_sig(iA, jj) = sum(A_alt_proj > 0)/length(A_alt_proj); % get the percentage of assemblies exceeding the activation cut off.
        J_sig(iA, jj) = sum(A_alt_proj > 0)/length(A_alt_proj); % get the percentage of assemblies exceeding the activation cut off.
        %
        % else
        %     % try a circ shift for sanity
        %     J_data = [];
        %     % for iS = size(A_out{jj}{1}.REM_pre_in, 2):-1:1
        %     %     J_data(:,iS) = circshift(A_out{jj}{1}.REM_post_in(:,iS), floor(MS_randn_range(1,1,1,length(J_data))));
        %     % end
        %     %
        %
        %     % for iShuff = 500:-1:1
        %     %     % temporary data from session jj;
        %     %     s_idx = randperm(size(A_out{jj}{1}.REM_post_in, 2));
        %     %     J_data = A_out{jj}{1}.REM_pre_in(:, s_idx);
        %     %
        %     %     % test if A_temps can be found in J-data
        %     %     T_proj = assembly_activity(A_temp,J_data');
        %     %
        %     %     J_proj = [];
        %     %     for tt = size(T_proj, 1):-1:1
        %     %         J_proj(tt) = sum(T_proj(tt,:)>9);
        %     %     end
        %     %     t_sigs(iShuff) = sum(J_proj > 0)/length(J_proj);
        %     %
        %     % end
        %     J_sig(iA, jj) = mean(t_sigs);
        % end
    end

    % t_data.Binary = A_out{iA}{1}.wake_data;
    % t_data.time = A_out{iA}{1}.wake_tvec;
    %
    %
    % [T_temp, T_proj, wake_data, wake_tvec, r_thresh] = MS_PCA_ICA_only(A_out{iA}{1}.wake_data, ones(1,length(A_out{iA}{1}.wake_data)), A_out{iA}{1}.bins,method, opts);
    %






end

%%
% try a circ shift for sanity
J_data = [];
% for iS = size(A_out{jj}{1}.REM_pre_in, 2):-1:1
%     J_data(:,iS) = circshift(A_out{jj}{1}.REM_post_in(:,iS), floor(MS_randn_range(1,1,1,length(J_data))));
% end
%

% for iShuff = 500:-1:1
%     % temporary data from session jj;
%     s_idx = randperm(size(A_out{jj}{1}.REM_post_in, 2));
%     J_data = A_out{jj}{1}.REM_pre_in(:, s_idx);
%
%     % test if A_temps can be found in J-data
%     T_proj = assembly_activity(A_temp,J_data');
%
%     J_proj = [];
%     for tt = size(T_proj, 1):-1:1
%         J_proj(tt) = sum(T_proj(tt,:)>9);
%     end
%     t_sigs(iShuff) = sum(J_proj > 0)/length(J_proj);
%
% end
J_sig(iA, jj) = mean(t_sigs);


%% same thing but using the same reactivation method as the REM data

c_ord = MS_linspecer(5);

rng(123, 'twister'); % for reproducibility.

opts = [];
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 500;
opts.threshold.permutations_percentile= 95;
opts.threshold.number_of_permutations= 500;


% init vars for number of assemblies per session and the shuff stats
a_num = NaN(size(A_out));
shuff_num = a_num;
wake_nA = []; wake_pre_nA = []; wake_post_nA = [];
J_sig = []; J_n = []; J_n_pval = []; J_n_ass = [];  J_r_ass = []; 
for iA = length(A_out):-1:1

    % current assembly wieghts;
    A_temp = A_out{iA}{1}.P_temp;
    % A_REM_post_proj = A_out{iA}{1}.REM_Post_proj;
    A_REM_tvec = A_out{iA}{1}.REM_Post_tvec; 


    % loop over sessions
    for jj = length(A_out):-1:1

        % J_data = ; % grab the data for the jjth session

        % test if A_temps can be found in J-data
        rng(123, 'twister'); % for reproducibility.
        A_alt_proj = assembly_activity(A_temp,A_out{jj}{1}.REM_Post_data');
        rng(123, 'twister'); % for reproducibility.
        A_REM_proj = assembly_activity(A_temp,A_out{iA}{1}.REM_Post_data');


        % trim the J_proj to be the same length as the A_proj
        if length(A_alt_proj) > length(A_REM_proj)
            A_alt_proj = A_alt_proj(:,1:length(A_REM_proj));
        end

        %get the same detection stats for both assemblies
        % % within session
        % rng(123, 'twister'); % for reproducibility.
        % [ReAct_stats, shuff.data, shuff.proj] = MS_Asmbly_proj_thresh(A_out{iA}{1}.REM_Post_data, A_temp, 500, 99);
        % ReAct_stats.p_val = [];
        % ReAct_stats.rate = [];
        % ReAct_stats.rate_p = [];
        % ReAct_stats.shuff_n = [];
        % ReAct_stats.shuff_r = [];
        % 
        % for ii = size(A_REM_post_proj,1):-1:1
        %     ReAct_stats.p_val(ii) = sum(sum(shuff.data > ReAct_stats.R_thresh,2) > sum(A_REM_post_proj(ii,:) > ReAct_stats.R_thresh))/ size(shuff.data,1);
        %     ReAct_stats.rate(ii) = sum(A_REM_post_proj(ii,:) > ReAct_stats.R_thresh) / ((A_out{iA}{1}.REM_Post_tvec(end) - A_out{iA}{1}.REM_Post_tvec(1))/60);
        %     ReAct_stats.shuff_rate = sum(shuff.data > ReAct_stats.R_thresh,2)./ ((A_out{iA}{1}.REM_Post_tvec(end) - A_out{iA}{1}.REM_Post_tvec(1))/60);
        %     ReAct_stats.rate_p(ii) = sum(ReAct_stats.shuff_rate > ReAct_stats.rate(ii)) / length(ReAct_stats.shuff_rate);
        % 
        %     ReAct_stats.shuff_n(ii) = median(sum(shuff.data > ReAct_stats.R_thresh,2));
        %     ReAct_stats.shuff_r(ii) = median(sum(shuff.data > ReAct_stats.R_thresh,2) / ((A_out{iA}{1}.REM_Post_tvec(end) - A_out{iA}{1}.REM_Post_tvec(1))/60));
        % 
        % end

        % alternative session
        rng(123, 'twister'); % for reproducibility.
        [Alt_stats, shuff.data, shuff.proj] = MS_Asmbly_proj_thresh(A_out{jj}{1}.REM_Post_data, A_temp, 500, 99);
        Alt_stats.p_val = [];
        Alt_stats.rate = [];
        Alt_stats.rate_p = [];
        Alt_stats.shuff_n = [];
        Alt_stats.shuff_r = [];

        for ii = size(A_alt_proj,1):-1:1
            Alt_stats.p_val(ii) = sum(sum(shuff.data > Alt_stats.R_thresh,2) > sum(A_alt_proj(ii,:) > Alt_stats.R_thresh))/ size(shuff.data,1);
            Alt_stats.rate(ii) = sum(A_alt_proj(ii,:) > Alt_stats.R_thresh) / ((A_REM_tvec(end) - A_REM_tvec(1))/60);
            Alt_stats.shuff_rate = sum(shuff.data > Alt_stats.R_thresh,2)./ ((A_REM_tvec(end) - A_REM_tvec(1))/60);
            Alt_stats.rate_p(ii) = sum(Alt_stats.shuff_rate > Alt_stats.rate(ii)) / length(Alt_stats.shuff_rate);

            Alt_stats.shuff_n(ii) = median(sum(shuff.data > Alt_stats.R_thresh,2));
            Alt_stats.shuff_r(ii) = median(sum(shuff.data > Alt_stats.R_thresh,2) / ((A_REM_tvec(end) - A_REM_tvec(1))/60));

        end

        % count the significant events per assembly on the jjth data
        % J_proj_cnt = []; A_proj_cnt = [];
        % for tt = size(A_alt_proj, 1):-1:1
        %     J_proj_cnt(tt) = nnz(A_alt_proj(tt,:)>ReAct_stats.R_thresh); % count the number of sig reactivations
        %     A_proj_cnt(tt) = nnz(A_REM_post_proj(tt,:)>ReAct_stats.R_thresh); % count the number of sig reactivations
        % end

        % J_sig(iA, jj) = sum(J_proj_cnt > 0)/length(J_proj_cnt); % get the percentage of assemblies exceeding the activation cut off.
        % J_n(iA, jj) = median(J_proj_cnt./A_proj_cnt); % get the percentage of assemblies exceeding the activation cut off.
        % J_n_pval(iA, jj) = sum(Alt_stats.p_val < 0.05)/sum(ReAct_stats.p_val < 0.05); % percent of jj assemblies passing the pval rest ./ number of real assemblies passing.
        J_n_ass(iA, jj) = sum(Alt_stats.p_val < 0.05); % number of jj assemblies passing the pval test
        J_r_ass(iA, jj) = median(Alt_stats.rate(Alt_stats.rate_p < 0.05)); % number of jj assemblies passing the pval test
        S_n_ass(iA, jj) = median(Alt_stats.shuff_n); % number of jj assemblies passing the pval test
        S_r_ass(iA, jj) = median(Alt_stats.shuff_r); % number of jj assemblies passing the pval test

        % J_n_ass(iA, jj,1) =sum(ReAct_stats.p_val < 0.05); % number of within assemblies passing the pval test
        %
            fprintf('iA: %d  | jj: %d\n', iA, jj)

    end
    sess_id{iA}  = [A_out{iA}{1}.info.subject '-' A_out{iA}{1}.info.session];
    
end
%% collect the number across conditions
idx = 1:size(J_n_ass,1); 
this_a_n = []; this_alt_n = []; 
this_a_r= []; this_alt_r = []; 
this_s_a_n= []; this_s_alt_n = []; 
this_s_a_r= []; this_s_alt_r = []; 

for ii = 1:size(J_n_ass,1)
    k_idx = idx ~=ii; 
    this_a_n(ii) = J_n_ass(ii,ii); 
    this_alt_n(ii) = median(J_n_ass(ii,k_idx), "omitnan"); 

    % get the rate metrics
    this_a_r(ii) = J_r_ass(ii,ii); 
    this_alt_r(ii) = median(J_r_ass(ii,k_idx),"omitnan"); 

    this_s_a_n(ii) = S_n_ass(ii,ii); 
    this_s_alt_n(ii) = median(S_n_ass(ii,k_idx), "omitnan"); 

    % get the rate metrics
    this_s_a_r(ii) = S_r_ass(ii,ii); 
    this_s_alt_r(ii) = median(S_r_ass(ii,k_idx), "omitnan"); 

end

%%
% quick stats of within vs alt values for each metric
d_idx = find(logical(eye(size(J_sig))));
off_idx = find(~logical(eye(size(J_sig))));

figure(104); clf
subplot(3,4,1)
imagesc(J_n_ass)
title('number of sig reactive assemblies')
cd = colorbar; 
cd.Position = cd.Position +[0.05 0 0 0];
c_a = clim; 
axis square
ylabel('session')
xlabel('session')

subplot(3,4,5)
[~,~,~,p, J_n_ass_stats] = MS_bar_w_err(this_a_n, this_alt_n,[c_ord(1,:); .7 .7 .7], 1, 'ttest');
set(gca, 'XTickLabel', {'Within' 'Across'})
ylabel('Num sig reactive assemblies')
y_n = ylim;

fprintf('Number of assemblies passing criteria within session (%0.2f +/- %0.2f) and across surrogate sessions (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.3f\n', ...
    mean(this_a_n), mean(this_alt_n), MS_SEM(this_a_n), MS_SEM(this_alt_n), ...
    J_n_ass_stats.df, J_n_ass_stats.tstat, p)

subplot(3,4,2)
imagesc(J_r_ass)
title('number of sig reactive assemblies')
cd = colorbar; 
cd.Position = cd.Position +[0.05 0 0 0];
c_r = clim; 
axis square


% same but for Rate
subplot(3,4,6)
[~,~,~,p, J_n_ass_stats] = MS_bar_w_err(this_a_r, this_alt_r,[c_ord(1,:); .7 .7 .7], 1, 'ttest');
set(gca, 'XTickLabel', {'Within' 'Across'})
title({'rate of reactive assemblies'})
y_r = ylim;

fprintf('Number of assemblies passing criteria within session (%0.2f +/- %0.2f) and across surrogate sessions (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.3f\n', ...
    mean(this_a_r), mean(this_alt_r), MS_SEM(this_a_r), MS_SEM(this_alt_r), ...
    J_n_ass_stats.df, J_n_ass_stats.tstat, p)

% SHUFFLE number
subplot(3,4,3)
imagesc(S_n_ass)
title({'Num sig reactive assemblies'; 'shuffle'})
cd = colorbar; 
cd.Position = cd.Position +[0.05 0 0 0];
clim(c_a)
axis square

subplot(3,4,7)
[~,~,~,p, S_n_ass_stats] = MS_bar_w_err(this_s_a_n, this_s_alt_n,[.3 .3 .3; .7 .7 .7], 1, 'ttest');
set(gca, 'XTickLabel', {'Within' 'Across'})
ylabel({'Num sig reactive assemblies'; 'shuffle'})
ylim(y_n);

fprintf('SHUFFLE Number of assemblies passing criteria within session (%0.2f +/- %0.2f) and across surrogate sessions (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.3f\n', ...
    mean(this_s_a_n), mean(this_s_alt_n), MS_SEM(this_s_a_n), MS_SEM(this_s_alt_n), ...
    S_n_ass_stats.df, S_n_ass_stats.tstat, p)

% SHUFFLE rate
subplot(3,4,4)
imagesc(S_r_ass)
title({'rate of reactive assemblies'; 'shuffle'})
cd = colorbar; 
clim(c_r)
cd.Position = cd.Position +[0.05 0 0 0];
axis square

subplot(3,4,8)
[~,~,~,p, S_r_ass_stats] = MS_bar_w_err(this_s_a_r, this_s_alt_r,[.3 .3 .3; .7 .7 .7], 1, 'ttest');
set(gca, 'XTickLabel', {'Within' 'Across'})
ylabel({'Reactivations / min'; 'shuffle'})
ylim(y_r);

fprintf('SHUFFLE Rate reactivations within session (%0.2f +/- %0.2f) and across surrogate sessions (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.3f\n', ...
    mean(this_s_a_r), mean(this_s_alt_r), MS_SEM(this_s_a_r), MS_SEM(this_s_alt_r), ...
    S_r_ass_stats.df, S_r_ass_stats.tstat, p)

% compare across vs shuffle

subplot(3,4,9)
[~,~,~,p, S_r_ass_stats] = MS_bar_w_err(this_a_n, this_s_a_n,[c_ord(1,:); .4 .4 .4], 1, 'ttest');
set(gca, 'XTickLabel', {'Within' 'Shuffle'})
title('number of sig reactive assemblies')
ylim(y_n);

fprintf('Num reactivations within session (%0.2f +/- %0.2f) and shuffles (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.3f\n', ...
    mean(this_a_n), mean(this_s_a_n), MS_SEM(this_a_n), MS_SEM(this_s_a_n), ...
    S_r_ass_stats.df, S_r_ass_stats.tstat, p)

subplot(3,4,10)
[~,~,~,p, S_r_ass_stats] = MS_bar_w_err(this_a_r, this_s_a_r,[c_ord(1,:); .4 .4 .4], 1, 'ttest');
set(gca, 'XTickLabel', {'Within' 'Shuffle'})
ylabel({'Reactivations / min'; 'shuffle'})
ylim(y_r);

fprintf('Rate reactivations across surrogate session (%0.2f +/- %0.2f) and shuffles (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.3f\n', ...
    mean(this_a_r), mean(this_s_a_r), MS_SEM(this_a_r), MS_SEM(this_s_a_r), ...
    S_r_ass_stats.df, S_r_ass_stats.tstat, p)


subplot(3,4,11)
[~,~,~,p, S_r_ass_stats] = MS_bar_w_err(this_alt_n, this_s_alt_n,[.7 .7 .7; .4 .4 .4], 1, 'ttest');
set(gca, 'XTickLabel', {'Across' 'Shuffle'})
title('number of sig reactive assemblies')
ylim(y_n);

fprintf('Num reactivations across surrogate session (%0.2f +/- %0.2f) and shuffles (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.3f\n', ...
    mean(this_alt_n), mean(this_s_alt_n), MS_SEM(this_alt_n), MS_SEM(this_s_alt_n), ...
    S_r_ass_stats.df, S_r_ass_stats.tstat, p)

subplot(3,4,12)
[~,~,~,p, S_r_ass_stats] = MS_bar_w_err(this_alt_r, this_s_alt_r,[.7 .7 .7; .4 .4 .4], 1, 'ttest');
set(gca, 'XTickLabel', {'Across' 'Shuffle'})
ylabel({'Reactivations / min'; 'shuffle'})
ylim(y_r);

fprintf('Rate reactivations across surrogate session (%0.2f +/- %0.2f) and shuffles (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.3f\n', ...
    mean(this_alt_r), mean(this_s_alt_r), MS_SEM(this_alt_r), MS_SEM(this_s_alt_r), ...
    S_r_ass_stats.df, S_r_ass_stats.tstat, p)

cfg_fig = []; 
cfg_fig.ft_size = 8; 
SetFigure(cfg_fig, gcf, 1); 
colormap(viridis)

%%
print("-bestfit",['C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\PC9\256_checks' filesep 'figS2_asmbly_checks_post'], '-dpdf', "-vector")

%%  PRE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(123, 'twister'); % for reproducibility.
opts = [];
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 500;
opts.threshold.permutations_percentile= 95;
opts.threshold.number_of_permutations= 500;

% init vars for number of assemblies per session and the shuff stats
a_num = NaN(size(A_out));
wake_nA = []; wake_pre_nA = []; wake_post_nA = [];
J_sig = []; J_n = []; J_n_pval = []; J_n_ass = [];  J_r_ass = []; 

for iA = length(A_out):-1:1

    % current assembly wieghts;
    A_temp = A_out{iA}{1}.P_temp;
    A_REM_tvec = A_out{iA}{1}.REM_Pre_tvec; 


    % loop over sessions
    for jj = length(A_out):-1:1

        % test if A_temps can be found in J-data
        rng(123, 'twister'); % for reproducibility.
        A_alt_proj = assembly_activity(A_temp,A_out{jj}{1}.REM_Pre_data');
        rng(123, 'twister'); % for reproducibility.
        A_REM_proj = assembly_activity(A_temp,A_out{iA}{1}.REM_Pre_data');

        % trim the J_proj to be the same length as the A_proj
        if length(A_alt_proj) > length(A_REM_proj)
            A_alt_proj = A_alt_proj(:,1:length(A_REM_proj));
        end

        % alternative session
        rng(123, 'twister'); % for reproducibility.
        [Alt_stats, shuff.data, shuff.proj] = MS_Asmbly_proj_thresh(A_out{jj}{1}.REM_Pre_data, A_temp, 500, 99);
        Alt_stats.p_val = [];
        Alt_stats.rate = [];
        Alt_stats.rate_p = [];
        Alt_stats.shuff_n = [];
        Alt_stats.shuff_r = [];

        for ii = size(A_alt_proj,1):-1:1
            Alt_stats.p_val(ii) = sum(sum(shuff.data > Alt_stats.R_thresh,2) > sum(A_alt_proj(ii,:) > Alt_stats.R_thresh))/ size(shuff.data,1);
            Alt_stats.rate(ii) = sum(A_alt_proj(ii,:) > Alt_stats.R_thresh) / ((A_REM_tvec(end) - A_REM_tvec(1))/60);
            Alt_stats.shuff_rate = sum(shuff.data > Alt_stats.R_thresh,2)./ ((A_REM_tvec(end) - A_REM_tvec(1))/60);
            Alt_stats.rate_p(ii) = sum(Alt_stats.shuff_rate > Alt_stats.rate(ii)) / length(Alt_stats.shuff_rate);

            Alt_stats.shuff_n(ii) = median(sum(shuff.data > Alt_stats.R_thresh,2));
            Alt_stats.shuff_r(ii) = median(sum(shuff.data > Alt_stats.R_thresh,2) / ((A_REM_tvec(end) - A_REM_tvec(1))/60));

        end

        J_n_ass(iA, jj) = sum(Alt_stats.p_val < 0.05); % number of jj assemblies passing the pval test
        J_r_ass(iA, jj) = median(Alt_stats.rate(Alt_stats.rate_p < 0.05)); % number of jj assemblies passing the pval test
        S_n_ass(iA, jj) = median(Alt_stats.shuff_n); % number of jj assemblies passing the pval test
        S_r_ass(iA, jj) = median(Alt_stats.shuff_r); % number of jj assemblies passing the pval test

            fprintf('iA: %d  | jj: %d\n', iA, jj)

    end
    sess_id{iA}  = [A_out{iA}{1}.info.subject '-' A_out{iA}{1}.info.session];
    
end
%% collect the number across conditions
idx = 1:size(J_n_ass,1); 
this_a_n = []; this_alt_n = []; 
this_a_r= []; this_alt_r = []; 
this_s_a_n= []; this_s_alt_n = []; 
this_s_a_r= []; this_s_alt_r = []; 

for ii = 1:size(J_n_ass,1)
    k_idx = idx ~=ii; 
    this_a_n(ii) = J_n_ass(ii,ii); 
    this_alt_n(ii) = median(J_n_ass(ii,k_idx), "omitnan"); 

    % get the rate metrics
    this_a_r(ii) = J_r_ass(ii,ii); 
    this_alt_r(ii) = median(J_r_ass(ii,k_idx),"omitnan"); 

    this_s_a_n(ii) = S_n_ass(ii,ii); 
    this_s_alt_n(ii) = median(S_n_ass(ii,k_idx), "omitnan"); 

    % get the rate metrics
    this_s_a_r(ii) = S_r_ass(ii,ii); 
    this_s_alt_r(ii) = median(S_r_ass(ii,k_idx), "omitnan"); 
end

%%
% quick stats of within vs alt values for each metric
d_idx = find(logical(eye(size(J_sig))));
off_idx = find(~logical(eye(size(J_sig))));

figure(104); clf
subplot(3,4,1)
imagesc(J_n_ass)
title('number of sig reactive assemblies')
cd = colorbar; 
cd.Position = cd.Position +[0.05 0 0 0];
c_a = clim; 
axis square
ylabel('session')
xlabel('session')

subplot(3,4,5)
[~,~,~,p, J_n_ass_stats] = MS_bar_w_err(this_a_n, this_alt_n,[c_ord(1,:); .7 .7 .7], 1, 'ttest');
set(gca, 'XTickLabel', {'Within' 'Across'})
ylabel('Num sig reactive assemblies')
y_n = ylim;

fprintf('Number of assemblies passing criteria within session (%0.2f +/- %0.2f) and across surrogate sessions (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.4f\n', ...
    mean(this_a_n), mean(this_alt_n), MS_SEM(this_a_n), MS_SEM(this_alt_n), ...
    J_n_ass_stats.df, J_n_ass_stats.tstat, p)

subplot(3,4,2)
imagesc(J_r_ass)
title('number of sig reactive assemblies')
cd = colorbar; 
cd.Position = cd.Position +[0.05 0 0 0];
c_r = clim; 
axis square


% same but for Rate
subplot(3,4,6)
[~,~,~,p, J_n_ass_stats] = MS_bar_w_err(this_a_r, this_alt_r,[c_ord(1,:); .7 .7 .7], 1, 'ttest');
set(gca, 'XTickLabel', {'Within' 'Across'})
title({'rate of reactive assemblies'})
y_r = ylim;

fprintf('Number of assemblies passing criteria within session (%0.2f +/- %0.2f) and across surrogate sessions (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.4f\n', ...
    mean(this_a_r), mean(this_alt_r), MS_SEM(this_a_r), MS_SEM(this_alt_r), ...
    J_n_ass_stats.df, J_n_ass_stats.tstat, p)

% SHUFFLE number
subplot(3,4,3)
imagesc(S_n_ass)
title({'Num sig reactive assemblies'; 'shuffle'})
cd = colorbar; 
cd.Position = cd.Position +[0.05 0 0 0];
clim(c_a)
axis square

subplot(3,4,7)
[~,~,~,p, S_n_ass_stats] = MS_bar_w_err(this_s_a_n, this_s_alt_n,[.3 .3 .3; .7 .7 .7], 1, 'ttest');
set(gca, 'XTickLabel', {'Within' 'Across'})
ylabel({'Num sig reactive assemblies'; 'shuffle'})
ylim(y_n);

fprintf('SHUFFLE Number of assemblies passing criteria within session (%0.2f +/- %0.2f) and across surrogate sessions (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.4f\n', ...
    mean(this_s_a_n), mean(this_s_alt_n), MS_SEM(this_s_a_n), MS_SEM(this_s_alt_n), ...
    S_n_ass_stats.df, S_n_ass_stats.tstat, p)

% SHUFFLE rate
subplot(3,4,4)
imagesc(S_r_ass)
title({'rate of reactive assemblies'; 'shuffle'})
cd = colorbar; 
clim(c_r)
cd.Position = cd.Position +[0.05 0 0 0];
axis square

subplot(3,4,8)
[~,~,~,p, S_r_ass_stats] = MS_bar_w_err(this_s_a_r, this_s_alt_r,[.3 .3 .3; .7 .7 .7], 1, 'ttest');
set(gca, 'XTickLabel', {'Within' 'Across'})
ylabel({'Reactivations / min'; 'shuffle'})
ylim(y_r);

fprintf('SHUFFLE Rate reactivations within session (%0.2f +/- %0.2f) and across surrogate sessions (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.4f\n', ...
    mean(this_s_a_r), mean(this_s_alt_r), MS_SEM(this_s_a_r), MS_SEM(this_s_alt_r), ...
    S_r_ass_stats.df, S_r_ass_stats.tstat, p)

% compare across vs shuffle

subplot(3,4,9)
[~,~,~,p, S_r_ass_stats] = MS_bar_w_err(this_a_n, this_s_a_n,[c_ord(1,:); .4 .4 .4], 1, 'ttest');
set(gca, 'XTickLabel', {'Within' 'Shuffle'})
title('number of sig reactive assemblies')
ylim(y_n);

fprintf('Num reactivations within session (%0.2f +/- %0.2f) and shuffles (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.4f\n', ...
    mean(this_a_n), mean(this_s_a_n), MS_SEM(this_a_n), MS_SEM(this_s_a_n), ...
    S_r_ass_stats.df, S_r_ass_stats.tstat, p)

subplot(3,4,10)
[~,~,~,p, S_r_ass_stats] = MS_bar_w_err(this_a_r, this_s_a_r,[c_ord(1,:); .4 .4 .4], 1, 'ttest');
set(gca, 'XTickLabel', {'Within' 'Shuffle'})
ylabel({'Reactivations / min'; 'shuffle'})
ylim(y_r);

fprintf('Rate reactivations across surrogate session (%0.2f +/- %0.2f) and shuffles (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.4f\n', ...
    mean(this_a_r), mean(this_s_a_r), MS_SEM(this_a_r), MS_SEM(this_s_a_r), ...
    S_r_ass_stats.df, S_r_ass_stats.tstat, p)


subplot(3,4,11)
[~,~,~,p, S_r_ass_stats] = MS_bar_w_err(this_alt_n, this_s_alt_n,[.7 .7 .7; .4 .4 .4], 1, 'ttest');
set(gca, 'XTickLabel', {'Across' 'Shuffle'})
title('number of sig reactive assemblies')
ylim(y_n);

fprintf('Num reactivations across surrogate session (%0.2f +/- %0.2f) and shuffles (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.4f\n', ...
    mean(this_alt_n), mean(this_s_alt_n), MS_SEM(this_alt_n), MS_SEM(this_s_alt_n), ...
    S_r_ass_stats.df, S_r_ass_stats.tstat, p)

subplot(3,4,12)
[~,~,~,p, S_r_ass_stats] = MS_bar_w_err(this_alt_r, this_s_alt_r,[.7 .7 .7; .4 .4 .4], 1, 'ttest');
set(gca, 'XTickLabel', {'Across' 'Shuffle'})
ylabel({'Reactivations / min'; 'shuffle'})
ylim(y_r);

fprintf('Rate reactivations across surrogate session (%0.2f +/- %0.2f) and shuffles (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.4f\n', ...
    mean(this_alt_r), mean(this_s_alt_r), MS_SEM(this_alt_r), MS_SEM(this_s_alt_r), ...
    S_r_ass_stats.df, S_r_ass_stats.tstat, p)

cfg_fig = []; 
cfg_fig.ft_size = 8; 
SetFigure(cfg_fig, gcf, 1); 
colormap(viridis)

%%
print("-bestfit",['C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\PC9\256_checks' filesep 'figS2_asmbly_checks_pre'], '-dpdf', "-vector")






%%
figure(103); clf
subplot(2,4,1)
% imagesc(J_sig*100)
% title('% of assemblies exceeding reactivation threshold of 9')
% colorbar
% axis square
% set(gca, 'YTick', 1:18, 'YTickLabel', sess_id)
% 
% subplot(2,4,2)
% imagesc(J_n*100)
% title('median n alt assemblies ./ n session assemblies')
% colorbar
% axis square

% subplot(2,4,3)
imagesc(J_n_pval*100)
title('percentage of sig reactive assemblies')
colorbar
axis square
% set(gca, 'YDir', 'normal')

% quick stats of within vs alt values for each metric
d_idx = find(logical(eye(size(J_sig))));
off_idx = find(~logical(eye(size(J_sig))));


subplot(2,4,5)
% [~,~,~,J_sig_stats] = MS_bar_w_err(J_sig(d_idx)*100, J_sig(off_idx)*100,[c_ord(1,:); .7 .7 .7], 1, 'ttest2');
% set(gca, 'XTickLabel', {'Within' 'Across'})
% ylabel('percentage of reactive assemblies')
% 
% 
% subplot(2,4,6)
% [~,~,~,J_n_stats] = MS_bar_w_err(J_n(d_idx)*100, J_n(off_idx)*100,[c_ord(1,:); .7 .7 .7], 1, 'ttest2');
% set(gca, 'XTickLabel', {'Within' 'Across'})
% ylabel('median alt:true reactivations')


% subplot(2,4,7)
[~,~,~,J_n_pval_stats] = MS_bar_w_err(J_n_pval(d_idx)*100, J_n_pval(off_idx)*100,[c_ord(1,:); .7 .7 .7], 1, 'ttest2');
set(gca, 'XTickLabel', {'Within' 'Across'})
ylabel('percentage of sig reactive assemblies')
y_lim = ylim;

subplot(2,4,6)
[d, b] = hist(J_n_pval(off_idx)*100, 50);
bar(b, d, 'FaceColor',  [.7 .7 .7]);
% xlim([0 ])
view(90,90)
set(gca, 'XDir', 'reverse')
xlim(y_lim)

subplot(2,4,7)
[~,~,~,p, J_n_ass_stats] = MS_bar_w_err(this_r_n, this_alt_n,[c_ord(1,:); .7 .7 .7], 1, 'ttest');
set(gca, 'XTickLabel', {'Within' 'Across'})
ylabel('Num sig reactive assemblies')
y_lim = ylim;

fprintf('Number of assemblies passing criteria within session (%0.2f +/- %0.2f) and across surrogate sessions (%0.2f +/- %0.2f); t(%d) = %0.2f, p = %0.3f\n', ...
    mean(this_r_n), mean(this_alt_n), MS_SEM(this_r_n), MS_SEM(this_alt_n), ...
    J_n_ass_stats.df, J_n_ass_stats.tstat, p)

subplot(2,4,8)
cla
hold on
[d, b] = hist(this_alt_n, 0:2.5:20);
area(b, d, 'FaceColor',  [.7 .7 .7]);

[d, b] = hist(this_r_n, 0:2.5:20);
area(b, d, 'FaceColor',  c_ord(1,:));
% xlim([0 ])
view(90,90)
set(gca, 'XDir', 'reverse')
xlim([y_lim])
%% quick R threshold generation

rng(123,'twister')
w_thresh = [];
for iA = 1
    data_h = A_out{iA}{1}.wake_data;

    nShuff = 100;
    wake_shuff_mat = [];

    Ass_shuff = NaN(1,nShuff);
    for iS = nShuff:-1:1
        tic
        shuff_data = NaN(size(data_h));
        for ic = 1:size(data_h,2)
            shuff_data(:,ic) = circshift(data_h(:,ic), floor(MS_randn_range(1,1,1,size(data_h,1))));
        end

        this_ass = assembly_patterns(shuff_data');
        if ~isempty(this_ass)
            S_prog = assembly_activity(this_ass,shuff_data');

            wake_shuff_mat(iS,:) =  S_prog(1,:);
            keep_idx(iS) = 1;
        else
            wake_shuff_mat(iS,:) = NaN(1,length(shuff_data));
            keep_idx(iS) = 0;
        end
        %     for ii = size(this_ass,2):-1:1

        if sum(max(this_ass) > 0.2) >0
            Ass_shuff(iS) = sum(max(this_ass) > 0.2);
        else
            Ass_shuff(iS) = 0;
        end
        %     end
        fprintf('Shuff # %.0f found %.0f assemblies and took %2.2f seconds\n', iS, size(this_ass,2), toc)
    end

    shuff_stats.shuff_n = Ass_shuff;
    shuff_stats.mean = mean(Ass_shuff);
    shuff_stats.sd = std(Ass_shuff);
    shuff_stats.p95 = prctile(Ass_shuff, 95, 'all');
    shuff_stats.p99 = prctile(Ass_shuff, 99, 'all');



    w_thresh(iA) = prctile(wake_shuff_mat(wake_shuff_mat >0), 99, 'all');

end

%% try the jaccard metric


c_ord = MS_linspecer(5);

for iS = 1:3

    this_A = A_out{iS}{1}.P_temp;
    this_data = A_out{iS}{1}.REM_Post_data;

    this_proj_alt = nan(length(this_data), 1);
    this_R_proj = A_out{iS}{1}.REM_Post_proj;

    % convert A_temp to only sig values
    for iA = 1:size(this_A,2)

        %%
        idx = zscore(this_A(:,iA)) <1 ;
        this_A_pos =  this_A(:,iA);
        this_A_pos(idx)  = 0;

        this_A_pos(~idx)  = 1;

        % use the classic but with the negative weights set to zero

        % this_A_pos
        this_proj_aa = assembly_activity(this_A(:,iA) ,this_data');
        this_proj_alt = assembly_activity(this_A(~idx,iA) ,this_data(:,~idx)');


        % TRY AN alternative distance metric
        for ii = 1:length(this_data)
            this_proj_alt2(ii) = pdist2(this_A_pos', this_data(ii,:), "cosine");
            this_proj_alt2(ii) = 1 - this_proj_alt2(ii);

            % jaccard requires removing the zeros first?
            % a_idx = sum(zscore(this_A_pos(:,1)) > 0, 2) > 0;
            %
            % jac_A = this_A_pos(a_idx);
            % jac_data = this_data(ii,a_idx);
            % this_proj_alt(ii) = pdist2(this_A_pos', this_data(ii,:), "jaccard");

        end

        % plot the alternative projection and the original.
        figure(iS*100+iA)
        clf
        subplot(6,4,[1 5 9 13 17 21])
        hold on
        stem(this_A(:,iA), 'color', [.8 .8 .8 .2])
        a_idx = sum(zscore(this_A(:,iA)) > 0, 2) > 0;
        stem(find(a_idx), this_A(find(a_idx),iA), 'color',winter(1), 'MarkerFaceColor', winter(1))
        view(90,90)

        ax(1) = subplot(6,4,2:4);
        plot(this_proj_aa(1,:), 'color', c_ord(1,:));
        title('standard all wieghts')

        ax(2) = subplot(6,4,6:8);
        plot(this_proj_alt', 'color',c_ord(2,:));
        title('standard pos only wieghts')

        ax(3) = subplot(6,4,10:12);
        plot(this_proj_alt2', 'color',c_ord(2,:));
        title('cosine distance')

        ax(4) = subplot(6,4,[14:16 22:24]);
        imagesc(this_data(:, a_idx)')

        linkaxes(ax, 'x')
        xlim([1 length(this_proj_alt)])


    end
end


MS_asmbly_quick_plot(this_A_pos, 1-this_proj_alt',this_data,1, 1)


%%
iA = 15;
jj = 16;

A_temp = A_out{iA}{1}.P_temp;

for ii = 1:size(A_temp,2)

    % y_max = max([A_out{iA}{1}.REM_Post_proj(ii,:), A_out{jj}{1}.REM_Post_proj(ii,:)]);

    figure(ii)
    clf
    A_data = A_out{iA}{1}.REM_Post_data;

    % test if A_temps can be found in J-data
    rng(123, 'twister'); % for reproducibility.
    A_alt_proj = assembly_activity(A_temp,A_data');

    % test if A_temps can be found in J-data
    rng(123, 'twister'); % for reproducibility.
    J_proj = assembly_activity(A_temp,A_out{jj}{1}.REM_Post_data');

    y_max = max([A_alt_proj(ii,:), J_proj(ii,:)]);

    MS_asmbly_quick_plot(A_out{iA}{1}.P_temp, A_alt_proj,A_data,ii )

    subplot(2,4,2:4);
    ylim([0 y_max])


    figure(ii+100)
    clf


    MS_asmbly_quick_plot(A_out{iA}{1}.P_temp, J_proj,A_out{jj}{1}.REM_Post_data,ii )
    ax(1) = subplot(2,4,2:4);
    ylim([0 y_max])


    % figure(ii+20)
    % clf
    % for kk = size()
    % MS_asmbly_quick_plot(A_out{iA}{1}.P_temp, A_out{jj}{1}.REM_Post_proj,A_out{jj}{1}.REM_Post_data,ii )
    % ax(1) = subplot(2,4,2:4);
    % ylim([0 y_max])


end

%% same thing but jst against shuffle

c_ord = MS_linspecer(5);

rng(123, 'twister'); % for reproducibility.

opts = [];
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 500;
opts.threshold.permutations_percentile= 95;
opts.threshold.number_of_permutations= 500;


% init vars for number of assemblies per session and the shuff stats
a_num = NaN(size(A_out));
shuff_num = a_num;
 S_n_pval = []; S_n_shuff_pval = []; 
for iA = length(A_out):-1:1

   % test if A_temps can be found in J-data
     
        rng(123, 'twister'); % for reproducibility.
        A_REM_post_proj = assembly_activity(A_temp,A_out{iA}{1}.REM_Post_data');


        Alt_temp_shuff_pval = [];
    for ii = size(A_REM_post_proj,1):-1:1
        ReAct_stats.p_val(ii) = sum(sum(shuff.data > ReAct_stats.R_thresh,2) > sum(A_REM_post_proj(ii,:) > ReAct_stats.R_thresh))/ size(shuff.data,1);
        ReAct_stats.rate(ii) = sum(A_REM_post_proj(ii,:) > ReAct_stats.R_thresh) / ((A_out{iA}{1}.REM_Post_tvec(end) - A_out{iA}{1}.REM_Post_tvec(1))/60);
        ReAct_stats.shuff_rate = sum(shuff.data > ReAct_stats.R_thresh,2)./ ((A_out{iA}{1}.REM_Post_tvec(end) - A_out{iA}{1}.REM_Post_tvec(1))/60);
        ReAct_stats.rate_p(ii) = sum(ReAct_stats.shuff_rate > ReAct_stats.rate(ii)) / length(ReAct_stats.shuff_rate);

        for jj = size(all_shuff_proj,2):-1:1
            Alt_temp_shuff_pval(ii,jj) = sum(sum(shuff.data > ReAct_stats.R_thresh,2) > sum(all_shuff_proj{jj}(ii,:) > ReAct_stats.R_thresh))/ size(shuff.data,1);
        end
        % Alt_stats.p_val(ii) = mean(Alt_temp_shuff_pval);
    end


  % % alternative session
  %       rng(123, 'twister'); % for reproducibility.
  %       [Alt_stats, shuff.data, shuff.proj] = MS_Asmbly_proj_thresh(A_out{jj}{1}.REM_Post_data, A_temp, 500, 99);
  %       Alt_stats.p_val = [];
  %       Alt_stats.rate = [];
  %       Alt_stats.rate_p = [];
  % 
  %       for ii = size(A_alt_proj,1):-1:1
  %           Alt_stats.p_val(ii) = sum(sum(shuff.data > ReAct_stats.R_thresh,2) > sum(A_alt_proj(ii,:) > ReAct_stats.R_thresh))/ size(shuff.data,1);
  %           Alt_stats.rate(ii) = sum(A_alt_proj(ii,:) > ReAct_stats.R_thresh) / ((A_out{iA}{1}.REM_Post_tvec(end) - A_out{iA}{1}.REM_Post_tvec(1))/60);
  %           Alt_stats.shuff_rate = sum(shuff.data > ReAct_stats.R_thresh,2)./ ((A_out{iA}{1}.REM_Post_tvec(end) - A_out{iA}{1}.REM_Post_tvec(1))/60);
  %           Alt_stats.rate_p(ii) = sum(Alt_stats.shuff_rate > Alt_stats.rate(ii)) / length(ReAct_stats.shuff_rate);
  %       end

    
    % % count the significant events per assembly on the jjth data
    % J_proj_cnt = []; A_proj_cnt = [];
    % for tt = size(A_alt_proj, 1):-1:1
    %     J_proj_cnt(tt) = nnz(A_alt_proj(tt,:)>9); % count the number of sig reactivations
    %     A_proj_cnt(tt) = nnz(A_REM_post_proj(tt,:)>9); % count the number of sig reactivations
    % end
    % 
    % J_sig(iA, jj) = sum(J_proj_cnt > 0)/length(J_proj_cnt); % get the percentage of assemblies exceeding the activation cut off.
    % J_n(iA, jj) = median(J_proj_cnt./A_proj_cnt); % get the percentage of assemblies exceeding the activation cut off.
    S_n_pval(iA) = sum(Alt_stats.p_val < 0.05)/sum(ReAct_stats.p_val < 0.05); % percent of jj assemblies passing the pval rest ./ number of real assemblies passing.
    
    for jj = size(Alt_temp_shuff_pval,2):-1:1
        this_shuff(jj) = sum(Alt_temp_shuff_pval(:,jj) < 0.05)/sum(ReAct_stats.p_val < 0.05); 
    end
    S_n_shuff_pval(iA, jj) = median(this_shuff); 

    %
    sess_id{iA}  = [A_out{iA}{1}.info.subject '-' A_out{iA}{1}.info.session];
    fprintf('iA: %d  | jj: %d\n', iA, jj)
end

% plot
figure(191)
clf
subplot(1,2,1)
[~,~,~,J_n_pval_stats] = MS_bar_w_err(S_n_pval*100, S_n_shuff_pval*100,[c_ord(1,:); .7 .7 .7], 1, 'ttest2');
set(gca, 'XTickLabel', {'Real' 'Shuff'})
ylabel('percentage of sig reactive assemblies')
y_lim = ylim;

subplot(2,4,8)
[d, b] = hist(J_n_pval(off_idx)*100, 50);
bar(b, d, 'FaceColor',  [.7 .7 .7]);
% xlim([0 ])
view(90,90)
set(gca, 'XDir', 'reverse')
xlim([y_lim])

% % shuff data
%     % test if A_temps can be found in J-data
%
%     %make a shuffle data set
%     shuff_data= [];
%     for ii = 10:-1:1
%         for jj = 1:size(A_out{iA}{1}.REM_Post_data,2)
%         shuff_data(:,jj,ii) = circshift(A_out{iA}{1}.REM_Post_data(:,jj), round(MS_randn_range(1, 1, 1, size(A_out{iA}{1}.REM_Post_data,1))));
%         end
%     end
%
%     % loop over the shuffles
%     for ii = 1:10
%         rng(123, 'twister'); % for reproducibility.
%         A_alt_proj(:,:,ii) = assembly_activity(A_temp,shuff_data(:,:,ii)');
%
%         % trim the J_proj to be the same length as the A_proj
%         if length(A_alt_proj) > length(A_REM_post_proj)
%             A_alt_proj = A_alt_proj(:,1:length(A_REM_post_proj),:);
%         end
%
%
%         % alternative session
%         rng(123, 'twister'); % for reproducibility.
%         [Alt_stats, shuff.data, shuff.proj] = MS_Asmbly_proj_thresh(A_out{jj}{1}.REM_Post_data, A_temp, 500, 99);
%         Alt_stats.p_val = [];
%         Alt_stats.rate = [];
%         Alt_stats.rate_p = [];
%
%         for ii = size(A_alt_proj,1):-1:1
%             Alt_stats.p_val(ii) = sum(sum(shuff.data > ReAct_stats.R_thresh,2) > sum(A_alt_proj(ii,:) > ReAct_stats.R_thresh))/ size(shuff.data,1);
%             Alt_stats.rate(ii) = sum(A_alt_proj(ii,:) > ReAct_stats.R_thresh) / ((A_out{iA}{1}.REM_Post_tvec(end) - A_out{iA}{1}.REM_Post_tvec(1))/60);
%             Alt_stats.shuff_rate = sum(shuff.data > ReAct_stats.R_thresh,2)./ ((A_out{iA}{1}.REM_Post_tvec(end) - A_out{iA}{1}.REM_Post_tvec(1))/60);
%             Alt_stats.rate_p(ii) = sum(Alt_stats.shuff_rate > Alt_stats.rate(ii)) / length(ReAct_stats.shuff_rate);
%         end
%
%     end
%


%% shuffle data and ICA

c_ord = MS_linspecer(5);

rng(123, 'twister'); % for reproducibility.

opts = [];
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 500;
opts.threshold.permutations_percentile= 95;
opts.threshold.number_of_permutations= 500;


% init vars for number of assemblies per session and the shuff stats
shuff_num = 500;

for iA = length(A_out):-1:1

    % grad some data
    A_temp = A_out{iA}{1}.P_temp;
    A_data = A_out{iA}{1}.wake_data; 

    A_REM_post_proj = A_out{iA}{1}.REM_Post_proj;



    % test if A_temps can be found in J-data
    rng(123, 'twister'); % for reproducibility.
    A_alt_proj = assembly_activity(A_temp,A_out{iA}{1}.REM_Post_data');


end