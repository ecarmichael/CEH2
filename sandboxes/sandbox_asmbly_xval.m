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

      

        if jj == iA
            J_data = A_out{iA}{1}.REM_post_in;

            % test if A_temps can be found in J-data
            T_proj = assembly_activity(A_temp,J_data');

            J_proj = [];
            for tt = size(T_proj, 1):-1:1
                J_proj(tt) = sum(T_proj(tt,:)>9);
            end

            J_sig(jj) = sum(J_proj > 0)/length(J_proj); % get the percentage of assemblies exceeding the activation cut off.

        else
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
        end
    end

    % t_data.Binary = A_out{iA}{1}.wake_data;
    % t_data.time = A_out{iA}{1}.wake_tvec;
    %
    %
    % [T_temp, T_proj, wake_data, wake_tvec, r_thresh] = MS_PCA_ICA_only(A_out{iA}{1}.wake_data, ones(1,length(A_out{iA}{1}.wake_data)), A_out{iA}{1}.bins,method, opts);
    %






end

%% same thing but using the same reactivation method as the REM data


for iA = length(A_out):-1:1




    for jj = length(A_out):-1:1




    end

end



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
        stem(this_A(:,1), 'color', [.8 .8 .8 .2])
        a_idx = sum(zscore(this_A(:,1)) > 0, 2) > 0;
        stem(find(a_idx), this_A(find(a_idx),1), 'color',winter(1), 'MarkerFaceColor', winter(1))
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

for ii = 1:5

    y_max = max([A_out{iA}{1}.REM_Post_proj(ii,:), A_out{jj}{1}.REM_Post_proj(ii,:)]); 

    figure(ii)
    clf
             J_data = A_out{iA}{1}.REM_post_in;

            % test if A_temps can be found in J-data
            T_proj = assembly_activity(A_temp,J_data');

    MS_asmbly_quick_plot(A_out{iA}{1}.P_temp, T_proj,J_data,ii )
    % MS_asmbly_quick_plot(A_out{iA}{1}.P_temp, A_out{iA}{1}.REM_Post_proj,A_out{iA}{1}.REM_Post_data,ii )
    subplot(2,4,2:4);
    ylim([0 y_max])

    figure(ii+10)
    clf
    J_data = A_out{jj}{1}.REM_post_in;

    % test if A_temps can be found in J-data
    T_proj = assembly_activity(A_temp,J_data');



    MS_asmbly_quick_plot(A_out{iA}{1}.P_temp, T_proj,J_data,ii )
    ax(1) = subplot(2,4,2:4);
    ylim([0 y_max])


    % figure(ii+20)
    % clf
    % for kk = size()
    % MS_asmbly_quick_plot(A_out{iA}{1}.P_temp, A_out{jj}{1}.REM_Post_proj,A_out{jj}{1}.REM_Post_data,ii )
    % ax(1) = subplot(2,4,2:4);
    % ylim([0 y_max])


end