%% sandbox_NAMI_SWR_analyses



%% loop over all the sessions of interest. Ultimately you will be able to loop over the sessions and hold the metrics of interest in order to get some session/subject level comparisons.
f_list = dir('pox*.mat');
k_idx = ones(length(f_list),1);


for iS = 1:length(f_list)
    % fname = 'pox3568_TFCD4'; % 3568 TFCD2/4/5 are all nice.


    % start with a nice session and generate some simple plots.

    this_sess = load([f_list(iS).name]);
    fname = fieldnames(this_sess.data);
    % Extract the session data for further analysis
    this_sess = this_sess.data.(fname{1}); % Access the first field of the session data

    % if length(strfind(fname{1}, '_')) > 1
    numText = regexp(fname{1}, '-?\d+\.?\d*', 'match', 'once');
    u_idx = strfind(fname{1}, '_'); 

        out.sess{iS} = upper(fname{1}(u_idx(end)+1:end));
        out.sub{iS} = ['pox' numText];

        if contains(out.sess{iS}, 'TL')
            out.sess{iS} = strrep(out.sess{iS}, 'TL', 'LT'); 
        end

    % end


    %%
    if isfield(this_sess, 'S')
        disp(fname)
    end

    disp(length(this_sess.swrs_ca1.tstart) )

    if length(this_sess.swrs_sub.tstart) < 200
        k_idx(iS) = 0;
        continue
    end
    %% Get the relative timing of each Sub SWR compared to the nearest Ca1 SWR
    % get the center of the events
    sub_cent = IVcenters(this_sess.swrs_sub);
    ca1_cent = IVcenters(this_sess.swrs_ca1);

    % loop over the events and find pairs and their offest

    pair_swr = []; % make an empty place to put things.
    offset_max = .75; % max time between CA1 peak and Sub peak + or minus.

    for ii = length(sub_cent):-1:1
        this_d = sub_cent(ii) - ca1_cent;
        co_evt = find(abs(this_d) < offset_max);

        if isempty(co_evt) || length(co_evt) > 1
            pair_swr(ii) = NaN;
        else
            pair_swr(ii) = this_d(co_evt);
        end
    end

    % check for why the number of std or atyp don't sum to total Ca1/Sub.

    % make a new set of events for Ca1 leading or Sub leading
    this_sess.swrs_sub_atyp = SelectIV([], this_sess.swrs_sub, pair_swr < -0.01);
    this_sess.swrs_sub_std = SelectIV([], this_sess.swrs_sub, pair_swr > 0.01);

    % find the CA1 SWRs that overlap with each Sub SWR type
    this_sess.swrs_ca1_atyp  = IntersectIV([], this_sess.swrs_ca1, this_sess.swrs_sub_atyp);
    this_sess.swrs_ca1_std  = IntersectIV([], this_sess.swrs_ca1, this_sess.swrs_sub_std);

    out.swrs_sub_atyp{iS} = this_sess.swrs_sub_atyp;
    out.swrs_sub_std{iS} = this_sess.swrs_sub_std;
    out.swrs_ca1_atyp{iS} = this_sess.swrs_ca1_atyp;
    out.swrs_ca1_std{iS} = this_sess.swrs_ca1_std;

    % basic measures
    % duration
    out.swr_sub_std_dur{iS} = out.swrs_sub_std{iS}.tend - out.swrs_sub_std{iS}.tstart;
    out.swr_sub_atyp_dur{iS} = out.swrs_sub_atyp{iS}.tend - out.swrs_sub_atyp{iS}.tstart;

    %


    %% quantify the participation of cells in the SWRs

    if isfield(this_sess, 'S')
        disp(fname)
        % std
        S_out_ca1_std =[];
        for ii = length(this_sess.swrs_ca1_std.tstart):-1:1
            this_S = restrict(this_sess.S, this_sess.swrs_ca1_std.tstart(ii), this_sess.swrs_ca1_std.tend(ii));
            dur = this_sess.swrs_ca1_std.tend(ii) - this_sess.swrs_ca1_std.tstart(ii);
            for jj = length(this_S.t):-1:1
                S_out_ca1_std(ii,jj) = length(this_S.t{jj}) / dur;
            end
        end

        % std
        S_out_sub_std =[];
        for ii = length(this_sess.swrs_sub_std.tstart):-1:1
            this_S = restrict(this_sess.S, this_sess.swrs_sub_std.tstart(ii), this_sess.swrs_sub_std.tend(ii));
            dur = this_sess.swrs_sub_std.tend(ii) - this_sess.swrs_sub_std.tstart(ii);
            for jj = length(this_S.t):-1:1
                S_out_sub_std(ii,jj) = length(this_S.t{jj}) / dur;
            end
        end

        % atyp
        S_out_sub_atyp =[];
        for ii = length(this_sess.swrs_sub_atyp.tstart):-1:1
            this_S = restrict(this_sess.S, this_sess.swrs_sub_atyp.tstart(ii), this_sess.swrs_sub_atyp.tend(ii));
            dur = this_sess.swrs_sub_atyp.tend(ii) - this_sess.swrs_sub_atyp.tstart(ii);
            for jj = length(this_S.t):-1:1
                S_out_sub_atyp(ii,jj) = length(this_S.t{jj}) /dur;
            end
        end

        % get the mean values
        out.Spk_ca1{iS} = mean(S_out_ca1_std, 1); % not used here.
        out.Spk_sub_std{iS} = mean(S_out_sub_std, 1);
        out.Spk_sub_atyp{iS} = mean(S_out_sub_atyp, 1);
        out.Spk_loc{iS} = this_sess.S.loc;
        cfg_rate = [];
        cfg_rate.DetectGaps = 1;
        this_sess.S = MS_spike_rates(this_sess.S, this_sess.csc.tvec, 2);
        for ii = length(this_sess.S.usr):-1:1
            out.Spk_rate{iS}(ii) = this_sess.S.usr{ii}.rate;
        end

    end


    %% Quantify the occurence of both types of SWRS within the task phases.
    % same thing you had done previsously using the restrict function, but this
    % time instead of printing the result you will want to hold the values in a
    % variable for later comparison.  YOu will need to add in the standard ad
    % atypical swrs too.

    % I would do this an array with some structure like this (rows are
    % conditions (pre, post, reward,...) and columns are the 4 types of SWR.

    swr_counts = NaN(7,4); % empty it to start.
    swr_counts_labels = {'Ca1', 'Sub', 'Std', 'Atyp'}; % labels for later.


    % for the LT sessions only get the pre (first 5mins as pre) and the
    % post (last 5mins)
    if contains(out.sess{iS}, 'LT')
        swr_ca1_pre = restrict(this_sess.swrs_ca1, this_sess.csc.tvec(1), this_sess.csc.tvec(1)+(5*60));
        swr_sub_pre = restrict(this_sess.swrs_sub, this_sess.csc.tvec(1),this_sess.csc.tvec(1)+(5*60));
        swr_sub_std_pre = restrict(this_sess.swrs_sub_std, this_sess.csc.tvec(1), this_sess.csc.tvec(1)+(5*60));
        swr_sub_atyp_pre = restrict(this_sess.swrs_sub_atyp, this_sess.csc.tvec(1), this_sess.csc.tvec(1)+(5*60));

        fprintf('Ca1 SWRs in Pre: <strong>%d</strong> \n', length(swr_ca1_pre.tstart))
        fprintf('Sub SWRs in Pre: <strong>%d</strong> \n', length(swr_sub_pre.tstart))
        fprintf('std Sub SWRs in Pre: <strong>%d</strong> \n', length(swr_sub_std_pre.tstart))
        fprintf('atyp Sub SWRs in Pre: <strong>%d</strong> \n', length(swr_sub_atyp_pre.tstart))

        % fill in the table. Coule be done in one line using a concatenation, but
        % this keeps it easy.
        swr_counts(1,1) = length(swr_ca1_pre.tstart);
        swr_counts(1,2) = length(swr_sub_pre.tstart);
        swr_counts(1,3) = length(swr_sub_std_pre.tstart);
        swr_counts(1,4) = length(swr_sub_atyp_pre.tstart);
        % add a label for the rows
        swr_counts_rows{1} = 'Pre';


        % NAMI TO DO:  same as above for the other measures.
        % restrict data to the post period.
        % ca1
        swr_ca1_post = restrict(this_sess.swrs_ca1,  this_sess.csc.tvec(end) - (5*60), this_sess.csc.tvec(end));
        %sub
        swr_sub_post = restrict(this_sess.swrs_sub, this_sess.csc.tvec(end) - (5*60), this_sess.csc.tvec(end));
        swr_sub_std_post = restrict(this_sess.swrs_sub_std, this_sess.csc.tvec(end) - (5*60), this_sess.csc.tvec(end));
        swr_sub_atyp_post = restrict(this_sess.swrs_sub_atyp, this_sess.csc.tvec(end) - (5*60), this_sess.csc.tvec(end));

        fprintf('Ca1 SWRs in Post: <strong>%d</strong> \n', length(swr_ca1_post.tstart))
        fprintf('Sub SWRs in Post: <strong>%d</strong> \n', length(swr_sub_post.tstart))
        fprintf('std Sub SWRs in Post: <strong>%d</strong> \n', length(swr_sub_std_post.tstart))
        fprintf('atyp Sub SWRs in Post: <strong>%d</strong> \n',length(swr_sub_atyp_post.tstart))

        % fill in the table.
        swr_counts(2,1) = length(swr_ca1_post.tstart);
        swr_counts(2,2) = length(swr_sub_post.tstart);
        swr_counts(2,3) = length(swr_sub_std_post.tstart);
        swr_counts(2,4) = length(swr_sub_atyp_post.tstart);
        % add a label for the rows
        swr_counts_rows{2} = 'Post';


    else

    % restrict data to the pre task baseline
    % ca1
    swr_ca1_pre = restrict(this_sess.swrs_ca1, this_sess.csc.tvec(1), this_sess.evts.t{2}(1));
    %sub
    swr_sub_pre = restrict(this_sess.swrs_sub, this_sess.csc.tvec(1), this_sess.evts.t{2}(1));
    swr_sub_std_pre = restrict(this_sess.swrs_sub_std, this_sess.csc.tvec(1), this_sess.evts.t{2}(1));
    swr_sub_atyp_pre = restrict(this_sess.swrs_sub_atyp, this_sess.csc.tvec(1), this_sess.evts.t{2}(1));


    fprintf('Ca1 SWRs in Pre: <strong>%d</strong> \n', length(swr_ca1_pre.tstart))
    fprintf('Sub SWRs in Pre: <strong>%d</strong> \n', length(swr_sub_pre.tstart))
    fprintf('std Sub SWRs in Pre: <strong>%d</strong> \n', length(swr_sub_std_pre.tstart))
    fprintf('atyp Sub SWRs in Pre: <strong>%d</strong> \n', length(swr_sub_atyp_pre.tstart))

    % fill in the table. Coule be done in one line using a concatenation, but
    % this keeps it easy.
    swr_counts(1,1) = length(swr_ca1_pre.tstart);
    swr_counts(1,2) = length(swr_sub_pre.tstart);
    swr_counts(1,3) = length(swr_sub_std_pre.tstart);
    swr_counts(1,4) = length(swr_sub_atyp_pre.tstart);
    % add a label for the rows
    swr_counts_rows{1} = 'Pre';


    % NAMI TO DO:  same as above for the other measures.
    % restrict data to the post period.
    % ca1
    swr_ca1_post = restrict(this_sess.swrs_ca1, this_sess.evts.t{2}(2), this_sess.csc.tvec(end));
    %sub
    swr_sub_post = restrict(this_sess.swrs_sub, this_sess.evts.t{2}(2), this_sess.csc.tvec(end));
    swr_sub_std_post = restrict(this_sess.swrs_sub_std, this_sess.evts.t{2}(2), this_sess.csc.tvec(end));
    swr_sub_atyp_post = restrict(this_sess.swrs_sub_atyp, this_sess.evts.t{2}(2), this_sess.csc.tvec(end));

    fprintf('Ca1 SWRs in Post: <strong>%d</strong> \n', length(swr_ca1_post.tstart))
    fprintf('Sub SWRs in Post: <strong>%d</strong> \n', length(swr_sub_post.tstart))
    fprintf('std Sub SWRs in Post: <strong>%d</strong> \n', length(swr_sub_std_post.tstart))
    fprintf('atyp Sub SWRs in Post: <strong>%d</strong> \n',length(swr_sub_atyp_post.tstart))

    % fill in the table.
    swr_counts(2,1) = length(swr_ca1_post.tstart);
    swr_counts(2,2) = length(swr_sub_post.tstart);
    swr_counts(2,3) = length(swr_sub_std_post.tstart);
    swr_counts(2,4) = length(swr_sub_atyp_post.tstart);
    % add a label for the rows
    swr_counts_rows{2} = 'Post';


    % get the baseline periods (30s before any tone) inside the TFC
    dur = 30;
    tone_t_on = sort([this_sess.evts.t{4}(1:2:end); this_sess.evts.t{5}(1:2:end)]);
    %ca1
    swr_ca1_base = restrict(this_sess.swrs_ca1, tone_t_on-30, tone_t_on);
    %sub
    swr_sub_base = restrict(this_sess.swrs_sub, tone_t_on-30, tone_t_on);
    swr_sub_std_base = restrict(this_sess.swrs_sub_std, tone_t_on-30, tone_t_on);
    swr_sub_atyp_base = restrict(this_sess.swrs_sub_atyp, tone_t_on-30, tone_t_on);

    fprintf('Ca1 SWRs in baseline: <strong>%d</strong> \n', length(swr_ca1_base.tstart))
    fprintf('Sub SWRs in baseline: <strong>%d</strong> \n', length(swr_sub_base.tstart))
    fprintf('std Sub SWRs in baseline: <strong>%d</strong> \n', length(swr_sub_std_base.tstart))
    fprintf('atyp Sub SWRs in baseline: <strong>%d</strong> \n', length(swr_sub_atyp_base.tstart))

    % fill in the table.
    swr_counts(3,1) = length(swr_ca1_base.tstart);
    swr_counts(3,2) = length(swr_sub_base.tstart);
    swr_counts(3,3) = length(swr_sub_std_base.tstart);
    swr_counts(3,4) = length(swr_sub_atyp_base.tstart);
    % add a label for the rows
    swr_counts_rows{3} = 'Baseline';


    % NAMI To Do:  same as above for the specific behaviour periods.
    % get the tone 1 periods
    tone_t_1 = sort([this_sess.evts.t{4}(1:2:end)]) ;

    swr_ca1_t1 = restrict(this_sess.swrs_ca1, tone_t_1, tone_t_1 +20);
    swr_sub_t1 = restrict(this_sess.swrs_sub, tone_t_1, tone_t_1 +20);
    swr_sub_std_t1 = restrict(this_sess.swrs_sub_std, tone_t_1, tone_t_1 +20);
    swr_sub_atyp_t1 = restrict(this_sess.swrs_sub_atyp, tone_t_1, tone_t_1 +20);

    fprintf('Ca1 SWRs in tone 1: <strong>%d</strong> \n', length(swr_ca1_t1.tstart))
    fprintf('Sub SWRs in tone 1: <strong>%d</strong> \n', length(swr_sub_t1.tstart))
    fprintf('std Sub SWRs in tone 1: <strong>%d</strong> \n', length(swr_sub_std_t1.tstart))
    fprintf('atyp Sub SWRs in tone 1: <strong>%d</strong> \n', length(swr_sub_atyp_t1.tstart))

    % fill in the table.
    swr_counts(4,1) = length(swr_ca1_t1.tstart);
    swr_counts(4,2) = length(swr_sub_t1.tstart);
    swr_counts(4,3) = length(swr_sub_std_t1.tstart);
    swr_counts(4,4) = length(swr_sub_atyp_t1.tstart);
    % add a label for the rows
    swr_counts_rows{4} = 'Tone 1';

    % get the tone 2 periods

    tone_t_2 = sort([this_sess.evts.t{5}(1:2:end)]);

    swr_ca1_t2 = restrict(this_sess.swrs_ca1, tone_t_2, tone_t_2 + 20);
    swr_sub_t2 = restrict(this_sess.swrs_sub, tone_t_2, tone_t_2 + 20);
    swr_sub_std_t2 = restrict(this_sess.swrs_sub_std, tone_t_2, tone_t_2 +20);
    swr_sub_atyp_t2 = restrict(this_sess.swrs_sub_atyp, tone_t_2, tone_t_2 +20);

    fprintf('Ca1 SWRs in tone 2: <strong>%d</strong> \n', length(swr_ca1_t2.tstart))
    fprintf('Sub SWRs in tone 2: <strong>%d</strong> \n', length(swr_sub_t2.tstart))
    fprintf('std Sub SWRs in tone 2: <strong>%d</strong> \n', length(swr_sub_std_t2.tstart))
    fprintf('atyp Sub SWRs in tone 2: <strong>%d</strong> \n', length(swr_sub_atyp_t2.tstart))

    % fill in the table.
    swr_counts(5,1) = length(swr_ca1_t2.tstart);
    swr_counts(5,2) = length(swr_sub_t2.tstart);
    swr_counts(5,3) = length(swr_sub_std_t2.tstart);
    swr_counts(5,4) = length(swr_sub_atyp_t2.tstart);
    % add a label for the rows
    swr_counts_rows{5} = 'Tone 2';

    % tone 1 trace period (20 sec after end of tone)

    tone_t_1end = sort([this_sess.evts.t{4}(2:2:end)]);

    swr_ca1_t1trace = restrict(this_sess.swrs_ca1, tone_t_1end, tone_t_1end + 20);
    swr_sub_t1trace = restrict(this_sess.swrs_sub, tone_t_1end, tone_t_1end + 20);
    swr_sub_std_t1trace = restrict(this_sess.swrs_sub_std, tone_t_1end, tone_t_1end + 20);
    swr_sub_atyp_t1trace = restrict(this_sess.swrs_sub_atyp, tone_t_1end, tone_t_1end + 20);

    fprintf('Ca1 SWRs in tone 1 trace: <strong>%d</strong> \n', length(swr_ca1_t1trace.tstart))
    fprintf('Sub SWRs in tone 1 trace: <strong>%d</strong> \n', length(swr_sub_t1trace.tstart))
    fprintf('std Sub SWRs in tone 1 trace: <strong>%d</strong> \n', length(swr_sub_std_t1trace.tstart))
    fprintf('atyp Sub SWRs in tone 1 trace: <strong>%d</strong> \n', length(swr_sub_atyp_t1trace.tstart))

    % fill in the table.
    swr_counts(6,1) = length(swr_ca1_t1trace.tstart);
    swr_counts(6,2) = length(swr_sub_t1trace.tstart);
    swr_counts(6,3) = length(swr_sub_std_t1trace.tstart);
    swr_counts(6,4) = length(swr_sub_atyp_t1trace.tstart);
    % add a label for the rows
    swr_counts_rows{6} = 'Tone 1 trace';

    % tone 2 trace period

    tone_t_2end = sort([this_sess.evts.t{5}(2:2:end)]);

    swr_ca1_t2trace = restrict(this_sess.swrs_ca1, tone_t_2end, tone_t_2end + 20);
    swr_sub_t2trace = restrict(this_sess.swrs_sub, tone_t_2end, tone_t_2end + 20);
    swr_sub_std_t2trace = restrict(this_sess.swrs_sub_std, tone_t_2end, tone_t_2end + 20);
    swr_sub_atyp_t2trace = restrict(this_sess.swrs_sub_atyp, tone_t_2end, tone_t_2end + 20);

    fprintf('Ca1 SWRs in tone 2 trace: <strong>%d</strong> \n', length(swr_ca1_t2trace.tstart))
    fprintf('Sub SWRs in tone 2 trace: <strong>%d</strong> \n', length(swr_sub_t2trace.tstart))
    fprintf('std Sub SWRs in tone 2 trace: <strong>%d</strong> \n', length(swr_sub_std_t2trace.tstart))
    fprintf('atyp Sub SWRs in tone 2 trace: <strong>%d</strong> \n', length(swr_sub_atyp_t2trace.tstart))

    % fill in the table.
    swr_counts(7,1) = length(swr_ca1_t1trace.tstart);
    swr_counts(7,2) = length(swr_sub_t1trace.tstart);
    swr_counts(7,3) = length(swr_sub_std_t1trace.tstart);
    swr_counts(7,4) = length(swr_sub_atyp_t1trace.tstart);
    % add a label for the rows
    swr_counts_rows{7} = 'Tone 2 trace';

    end
    % convert to proportions

    out.swr_prop{iS} = swr_counts ./ (swr_counts(:,3) + swr_counts(:,4));
    out.swr_rows{iS} = swr_counts_rows;

    clearvars -except iS out f_list k_idx
end % end cross sessions.




%% intra session stats



% compile spike activity across atyp and std SWRs
Spk_std = cell2mat(out.Spk_sub_std);
Spk_atyp = cell2mat(out.Spk_sub_atyp);
Spk_loc = logical(cell2mat(out.Spk_loc));
Spk_pyr = cell2mat(out.Spk_rate) < 20;


% collect the behaviours.
hab_idx = contains(out.sess, 'TFCD1') | contains(out.sess, 'TFCD2');
train_idx = contains(out.sess, 'TFCD3');
test_idx = contains(out.sess, 'TFCD4') | contains(out.sess, 'TFCD5');
LT_idx = contains(out.sess, 'LT'); 
pox_idx = contains(out.sub, '3265'); 
k_idx = logical(k_idx'); 


% convert to a matrix
prop_mat = []; 
for ii = length(out.swr_prop):-1:1
    if ~k_idx(ii) || pox_idx(ii)
    prop_mat(:,:,ii) = NaN(7,4); 

    else
    prop_mat(:,:,ii) = out.swr_prop{ii}; 
    end
end


swr_prop.hab = prop_mat(:,:,~pox_idx & hab_idx & k_idx)*100; 
swr_prop.train = prop_mat(:,:,~pox_idx & train_idx & k_idx)*100; 
swr_prop.test = prop_mat(:,:,~pox_idx & test_idx & k_idx)*100; 
swr_prop.all = prop_mat(:,:,~pox_idx & ~LT_idx & k_idx)*100; 




%% Initialize some things: 


c_ord = MS_linspecer(5); %generate a few nice colours.

fig_size = [100, 100, 500, 400]; % for consistent figure sizes.
ft_size = 12;  %font size
plot_flag =1;

if ispc
    fig_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\Inter_data\HF_Sub_SWR\Nami_figs';
else
    fig_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Inter_data/HF_Sub_SWR/Nami_figs';
end

%% plot all of the data from across the subjects  %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if plot_flag == 1
    figure(301)
    clf
    subplot(2,2,1)
    MS_bar_w_err(Spk_std(Spk_loc& Spk_pyr),  Spk_atyp(Spk_loc& Spk_pyr), c_ord([1,3],:), 1, 'ttest', 1:2, 1);
    ylim([0 100])
    set(gca, 'XTickLabel', {'Std', 'Atyp'}, 'yscale', 'log','ytick', [0 1 10 100], 'YTickLabel', {'0','1', '10', '100'},'fontsize', 10)
    title('Ca1 Pyramidal Cells')
    ylabel({'Within Ripple';'Firing Rate (Hz)'})


    subplot(2,2,2)
    [hb, eb, sc, p, stats]= MS_bar_w_err(Spk_std(~Spk_loc& Spk_pyr),  Spk_atyp(~Spk_loc& Spk_pyr), c_ord([1,3],:), 1, 'ttest', 1:2, 1);
    ylim([0 100])
    set(gca, 'XTickLabel', {'Std', 'Atyp'}, 'yscale', 'log','ytick', [0 1 10 100], 'YTickLabel', {'0','1', '10', '100'},'fontsize', 10)
    title('Sub Pyramidal Cells')


    %
     subplot(2,2,3)
     MS_bar_w_err(Spk_std(Spk_loc & ~Spk_pyr),  Spk_atyp(Spk_loc & ~Spk_pyr), c_ord([1,3],:), 1, 'ttest', 1:2);
     ylim([0 200])
     set(gca, 'XTickLabel', {'Std', 'Atyp'}, 'yscale', 'log','ytick', [0 1 10 100], 'YTickLabel', {'0','1', '10', '100'},'fontsize', 10)
     title('Ca1 interneurons')
     ylabel({'Within Ripple';'Firing Rate (Hz)'})

     subplot(2,2,4)
    [hb, eb, sc, p, stats]= MS_bar_w_err(Spk_std(~Spk_loc& ~Spk_pyr),  Spk_atyp(~Spk_loc & ~Spk_pyr), c_ord([1,3],:), 1, 'ttest', 1:2);
    % hb.E
    hb.FaceColor = 'none';
        ylim([0 200])
     set(gca, 'XTickLabel', {'Std', 'Atyp'}, 'yscale', 'log','ytick', [0 1 10 100], 'YTickLabel', {'0','1', '10', '100'},'fontsize', 10)
     title('Sub interneurons')


    exportgraphics(gcf, [fig_dir filesep 'Spikes_in_SWRs.pdf'], 'ContentType', 'vector');
end


%%  plot the duration of the events

if plot_flag == 1

    figure(3004)
    clf
    subplot(221)

    MS_rain_plot([cell2mat(out.swr_sub_std_dur')' cell2mat(out.swr_sub_atyp_dur')']'*1000,[zeros(size(cell2mat(out.swr_sub_std_dur')))', ones(size(cell2mat(out.swr_sub_atyp_dur')))']',...
        c_ord([1,3],:), 'ttest2', 1:2);
    xlim([20 140]); ylim([0.5 2.75])
    set(gca, 'yTickLabel', {'Std', 'Atyp'}, 'xscale', 'linear','ytick', [1 2], 'YTickLabel', {'Std', 'Atyp'},'fontsize', 10)
    xlabel('Swr Duration (ms)')

    exportgraphics(gcf, [fig_dir filesep 'Sub_Swr_duration.pdf'], 'ContentType', 'vector');


end

%% get the proportion of SWRs per task phase

figure(3005)
clf

cond = 'all';
subplot(2,3,1)
cla
h1 = MS_bar_w_err(squeeze(mean(swr_prop.(cond)(1:2,3,:), 'omitmissing'))', squeeze(mean(swr_prop.(cond)(1:2,4,:), 'omitmissing'))', [c_ord(1,:); c_ord(1,:)], 1, 'ttest', 1:2, 1);
ylim([0 100])
set(gca, 'xtick', [1.5], 'XTickLabel', {'pre+post'})

legend(h1, 'Std', 'box', 'off')
title(['Sub SWRs in ' cond ' '])
ylabel('% of Sub SWRs')

cond = 'all';
subplot(2,3,4)
cla
h1 = MS_bar_w_err(squeeze(swr_prop.(cond)(1,4,:))', squeeze(swr_prop.(cond)(2,4,:))', [c_ord(1,:)*.75; c_ord(1,:)], 1, 'ttest', 1:2);
h1.EdgeColor = 'none';
ylim([0 100])
set(gca, 'xtick', [1 2], 'XTickLabel', {'pre', 'post'})
ylabel('% of Atypical Sub SWRs')
title(['Sub SWRs in ' cond ' '])


% overall
cond = 'all';
subplot(2,3,2)
cla
data1 = squeeze(mean(swr_prop.(cond)(1:2,4,:), 'omitmissing'))'; 
data2 = squeeze(swr_prop.(cond)(3,4,:))';
data3 = squeeze(mean(swr_prop.(cond)(4:5,4,:), 'omitmissing'))'; 
data4 = squeeze(mean(swr_prop.(cond)(6:7,4,:), 'omitmissing'))';

n_idx = isnan(sum([data1; data2; data3; data4])); 

MS_bar_w_err4(data1(~n_idx)',data2(~n_idx)', data3(~n_idx)', data4(~n_idx)',...
    c_ord(1:4,:), 1, 'anova2', 1:4, {'pre+post', 'baseline', 'tones', 'trace'});

ylim([0 100])
set(gca,'XTickLabel', {'pre+post', 'baseline', 'tones', 'trace'})

title(['Sub SWRs in ' cond ' '])
ylabel('% of Atypical Sub SWRs')

cond = 'hab';
subplot(2,3,3)
cla
data1 = squeeze(mean(swr_prop.(cond)(1:2,4,:), 'omitmissing'))'; 
data2 = squeeze(swr_prop.(cond)(3,4,:))';
data3 = squeeze(mean(swr_prop.(cond)(4:5,4,:), 'omitmissing'))'; 
data4 = squeeze(mean(swr_prop.(cond)(6:7,4,:), 'omitmissing'))';

n_idx = isnan(sum([data1; data2; data3; data4])); 

MS_bar_w_err4(data1(~n_idx)',data2(~n_idx)', data3(~n_idx)', data4(~n_idx)',...
    c_ord(1:4,:), 1, 'anova2', 1:4, {'pre+post', 'baseline', 'tones', 'trace'});

ylim([0 100])
set(gca,'XTickLabel', {'pre+post', 'baseline', 'tones', 'trace'})
title(['Sub SWRs in ' cond ' '])


cond = 'train';
subplot(2,3,5)
cla
MS_bar_w_err4(squeeze(mean(swr_prop.(cond)(1:2,4,:), 'omitmissing'))', squeeze(mean(swr_prop.(cond)(3,4,:), 'omitmissing'))',...
    squeeze(mean(swr_prop.(cond)(4:5,4,:), 'omitmissing'))', squeeze(mean(swr_prop.(cond)(6:7,4,:), 'omitmissing'))',...
    c_ord(1:4,:), 1, [], 1:4, {'pre+post', 'baseline', 'tones', 'trace'});

ylim([0 100])
set(gca,'XTickLabel', {'pre+post', 'baseline', 'tones', 'trace'})
title(['Sub SWRs in ' cond ' '])

cond = 'test';
subplot(2,3,6)
cla
data1 = squeeze(mean(swr_prop.(cond)(1:2,4,:), 'omitmissing'))'; 
data2 = squeeze(swr_prop.(cond)(3,4,:))';
data3 = squeeze(mean(swr_prop.(cond)(4:5,4,:), 'omitmissing'))'; 
data4 = squeeze(mean(swr_prop.(cond)(6:7,4,:), 'omitmissing'))';

n_idx = isnan(sum([data1; data2; data3; data4])); 

MS_bar_w_err4(data1(~n_idx)',data2(~n_idx)', data3(~n_idx)', data4(~n_idx)',...
    c_ord(1:4,:), 1, 'anova2', 1:4, {'pre+post', 'baseline', 'tones', 'trace'});

ylim([0 100])
set(gca,'XTickLabel', {'pre+post', 'baseline', 'tones', 'trace'})

title(['Sub SWRs in ' cond ' '])

    exportgraphics(gcf, [fig_dir filesep 'task_phase.pdf'], 'ContentType', 'vector');

%% make an event triggered spectrogram of the CA1 and Sub SWRs relative to the CA1 center using fieldtrip.

% skip this if you don't want to deal with fieldtrip.

% get things ready for field trip.  You will need to clone the fieldtrip
% repo and initialize it by moving to the fieldtrip folder in matlab and
% running ft_defaults.

% pull Ca1 channel out first.
if plot_flag
    csc_ft = this_sess.csc;
    csc_ft.data = this_sess.csc.data(1,:);
    csc_ft.label = [];
    csc_ft.label{1} = this_sess.csc.label{1};

    [csc_ft_out, TFR] = Triggered_Spec_FT(csc_ft, ca1_cent, [], 80:.5:200, [], [-.25 .25], 1);
    set(gcf,'Units','pixels','position',fig_size); % helpful to keep the figures the same size.

    xlim([-.05 .05]); ylabel('frequency (hz)'); xlabel('time from CA1 SWR center (ms)')
    set(gca,'xtick', -.05:0.025:.05);
    set(gca, 'XTickLabel', get(gca, 'XTick')*1000)
    title('Ca1 SWRS')
    axis("square")

    % move this to a subplot.
    % ax = gca;
    % figure(3000)
    % subplot(3,2,1, ax)

    % Next Sub aligned to CA1

    csc_ft = this_sess.csc;
    csc_ft.data = this_sess.csc.data(2,:);
    csc_ft.label = [];
    csc_ft.label{1} = this_sess.csc.label{2};

    [csc_ft_out, TFR] = Triggered_Spec_FT(csc_ft, ca1_cent, [], 80:.5:200, [], [-.25 .25], 1);
    set(gcf,'Units','pixels','position',fig_size); % helpful to keep the figures the same size.

    xlim([-.05 .05]); ylabel('frequency (hz)'); xlabel('time from CA1 SWR center (ms)')
    set(gca,'xtick', -.05:0.025:.05);
    set(gca, 'XTickLabel', get(gca, 'XTick')*1000)
    title('Sub SWRS')
    axis("square")


    % Next Sub LEF aligned to CA1 standard SWR
    ca1_std_cent = IVcenters(this_sess.swrs_ca1_std);

    [csc_ft_out, TFR] = Triggered_Spec_FT(csc_ft, ca1_std_cent, [], 80:.5:200, [], [-.25 .25], 1);
    set(gcf,'Units','pixels','position',fig_size); % helpful to keep the figures the same size.

    xlim([-.05 .05]); ylabel('frequency (hz)'); xlabel('time from CA1 SWR center (ms)')
    set(gca,'xtick', -.05:0.025:.05);
    set(gca, 'XTickLabel', get(gca, 'XTick')*1000)
    title('Sub SWRS - Standard')
    axis("square");


    % Next Sub LEF aligned to CA1 atypical SWR
    ca1_atyp_cent = IVcenters(this_sess.swrs_ca1_atyp);

    [csc_ft_out, TFR] = Triggered_Spec_FT(csc_ft, ca1_atyp_cent, [], 80:.5:200, [], [-.25 .25], 1);
    set(gcf,'Units','pixels','position',fig_size); % helpful to keep the figures the same size.

    xlim([-.05 .05]); ylabel('frequency (hz)'); xlabel('time from CA1 SWR center (ms)')
    set(gca,'xtick', -.05:0.025:.05);
    set(gca, 'XTickLabel', get(gca, 'XTick')*1000)
    title('Sub SWRS - Atypical')
    axis("square")



    %% make a figure showing the distribution of CA1-Sub SWR timing within a session.
    % subplot(3,2,5) % subplot
    figure(3002); clf
    set(gcf,'Units','pixels','position',fig_size); % helpful to keep the figures the same size.

    hold on
    all_pairs = pair_swr(~isnan(pair_swr))*1000; %take the non-nan pair times and converting them to miliseconds

    pos_pairs = all_pairs(all_pairs> 0);
    neg_pairs = all_pairs(all_pairs< 0);

    % make two histograms (one for positive and one for negative so that they
    % can have different colours

    histogram(pos_pairs, -25:1:25, 'FaceColor',c_ord(1,:), 'EdgeColor',c_ord(1,:)) % make a histogram  within +/- 25ms
    histogram(neg_pairs, -25:1:25, 'FaceColor',c_ord(2,:), 'EdgeColor',c_ord(2,:)) % make a histogram  within +/- 25ms

    xline(0, '--k', 'LineWidth',1)
    % axis('square') % make things square.

    % Calculate and display the mean and standard deviation of the offsets
    meanOffset = mean(pair_swr, 'omitmissing');
    stdOffset = std(pair_swr, 'omitmissing');
    title(sprintf('Mean: %.2f ms, Std: %.2f ms', meanOffset * 1000, stdOffset * 1000), 'FontSize',ft_size, 'FontWeight','normal');

    % set the font to be a consistent size.

    set(gca, 'FontSize',ft_size)
    axis("square")


    %% NAMI TODO: plot some examples of the filtering and SWR detection.
    % this will compliment the above plot with one subplot for CA1 and Sub each
    % with a hold on that keeps the filtered trace, the zscored amplitude, and
    % the 2sd threshold.
    %
    % you should have the elements you need in that filtering script you made a
    % while back.  In case you don't have access to that. Remember that the
    % abslute of the Hilbert transform [abs(hilbert(csc_f.data(1,:)))] is how
    % you extract the amplitude.
    %
    % bonus points for having the filtered traces match the colours of the MUA
    % in the above plot. Also making the amplitue plot thicker with a linewidth
    % of 2, and having the threshold be a dashed line.

    % filter into the SWR band csc_f.data(1,:) is CA1 and csc_f.data(2,:) will
    % be the subiculum channel.
    cfg = [];
    cfg.f = [125 200];
    csc_f = FilterLFP(cfg, csc);


    % keep the plot sizes consistent.
    figure(202)
    clf
    set(gcf,'Units','inch','OuterPosition',[1 6 4 8]);

    % plot the Ca1 data
    subplot(2,1,1)
    hold on

    xlim([181.5 182.5]) % x limits for the nice SWR in the above plot

    % plot the Subiculum data
    subplot(2,1,2)
    hold on

    xlim([181.5 182.5])

end


