%% sandbox_NAMI_MUA_analyses



%% loop over all the sessions of interest. Ultimately you will be able to loop over the sessions and hold the metrics of interest in order to get some session/subject level comparisons.
f_list = dir('pox*.mat');

S_idx =  [3 12 13 16 18]; 

for iS = 1:length(S_idx)
    % fname = 'pox3568_TFCD4'; % 3568 TFCD2/4/5 are all nice.


    % start with a nice session and generate some simple plots.

    this_sess = load([f_list(S_idx(iS)).name]);
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


    %% Get the relative timing of each Sub MUA compared to the nearest Ca1 SWR

    cfg_rate = [];
    cfg_rate.DetectGaps = 1;
    this_sess.S = MS_spike_rates(this_sess.S, this_sess.csc.tvec, 2);
    for ii = length(this_sess.S.usr):-1:1
        out.Spk_rate{iS}(ii) = this_sess.S.usr{ii}.rate;
    end

    pyr_idx = out.Spk_rate{iS} < 20; 

    cfg_mua = [];
    cfg_mua.tvec = this_sess.csc.tvec;
    cfg_mua.sigma = 0.005;
    
    % Ca1 MUA
    S_ca = this_sess.S;
    S_ca.t(~this_sess.S.loc | ~pyr_idx) = [];
    S_ca.label(~this_sess.S.loc| ~pyr_idx) = [];
    S_ca.usr(~this_sess.S.loc| ~pyr_idx) = [];
    mua = getMUA(cfg_mua, S_ca); % get the smoothed multiunit activity;
    % plot(mua.tvec, zscore(mua.data),'color', this_sess.S.c_ord(2,:))
    % convert to TS for peaks
    [~, pk_idx] = findpeaks(zscore(mua.data),mua.tvec,  'MinPeakDistance', .05, 'MinPeakHeight',4);
    ca1_ts = ts({pk_idx}, {'Ca1'}); 

    % Sub MUA
    S_sub = this_sess.S;
    S_sub.t(this_sess.S.loc | ~pyr_idx) = [];
    S_sub.label(this_sess.S.loc| ~pyr_idx) = [];
    S_sub.usr(this_sess.S.loc| ~pyr_idx) = [];
    mua_sub = getMUA(cfg_mua, S_sub); % get the smoothed multiunit activity;
        % plot(mua_sub.tvec, zscore(mua_sub.data),'color', this_sess.S.c_ord(2,:))
    % convert to TS for peaks
    [~, pk_idx] = findpeaks(zscore(mua_sub.data),mua_sub.tvec,  'MinPeakDistance', .05, 'MinPeakHeight',4);
    sub_ts = ts({pk_idx}, {'Sub'}); 


    % get the xcorr between Ca1 and Sub MUA bursts;
    cfg_ccf = []; 
    cfg_ccf.binsize = 0.001; 
    cfg_ccf.max_t = .05;
    cfg_ccf.smooth =1; 
    cfg_ccf.gauss_w = .025; % width of Gaussian convolution window (in s)
    cfg_ccf.gauss_sd = 0.005;
    [out.mua_ccf{iS}, out.t_out{iS}]  = ccf(cfg_ccf, ca1_ts.t{1}, sub_ts.t{1}); 

    pair_mua = []; c_len = []; 
    % split the TS based on std or atyp; 
    for ii = length(sub_ts.t{1}):-1:1
        this_d = sub_ts.t{1}(ii) - ca1_ts.t{1};
        co_evt = find(abs(this_d) < .05);
                c_len(ii) = length(co_evt); 

        if isempty(co_evt) || length(co_evt) > 1
            pair_mua(ii) = NaN;
        else
            pair_mua(ii) = this_d(co_evt);
        end
    end

    % split them out. 
    std_ts = sub_ts; 
    std_ts.t{1} = sub_ts.t{1}(pair_mua>=0);

    atyp_ts = sub_ts; 
    atyp_ts.t{1} = sub_ts.t{1}(pair_mua<0);

    out.ca1_ts{iS} = ca1_ts;
    out.sub_ts{iS} = sub_ts;
    out.std_ts{iS} = std_ts;
    out.atyp_ts{iS} = atyp_ts;



    %% quantify the participation of cells in the SWRs
    mua_std_iv = iv(std_ts.t{1}-.05, std_ts.t{1}-.05);
    mua_atyp_iv = iv(atyp_ts.t{1}-.05, atyp_ts.t{1}-.05);

        % std
        S_out_sub_std =[];
        for ii = length(mua_std_iv.tstart):-1:1
            this_S = restrict(this_sess.S, mua_std_iv.tstart(ii), mua_std_iv.tend(ii));
            dur = mua_std_iv.tend(ii) - mua_std_iv.tstart(ii);
            for jj = length(this_S.t):-1:1
                S_out_sub_std(ii,jj) = length(this_S.t{jj}) / dur;
            end
        end

        % atyp
        S_out_sub_atyp =[];
        for ii = length(mua_atyp_iv.tstart):-1:1
            this_S = restrict(this_sess.S, mua_atyp_iv.tstart(ii), mua_atyp_iv.tend(ii));
            dur = mua_atyp_iv.tend(ii) - mua_atyp_iv.tstart(ii);
            for jj = length(this_S.t):-1:1
                S_out_sub_atyp(ii,jj) = length(this_S.t{jj}) /dur;
            end
        end

        % get the mean values
        out.Spk_sub_std{iS} = mean(S_out_sub_std, 1);
        out.Spk_sub_atyp{iS} = mean(S_out_sub_atyp, 1);
        out.Spk_loc{iS} = this_sess.S.loc;

    %% Quantify the occurence of both types of SWRS within the task phases.
    % same thing you had done previsously using the restrict function, but this
    % time instead of printing the result you will want to hold the values in a
    % variable for later comparison.  YOu will need to add in the standard ad
    % atypical swrs too.

    % I would do this an array with some structure like this (rows are
    % conditions (pre, post, reward,...) and columns are the 4 types of SWR.

    mua_counts = NaN(7,4); % empty it to start.
    mua_counts_labels = {'Ca1', 'Sub', 'Std', 'Atyp'}; % labels for later.

    phase_labels = {'pre', 'post', 'baseline', 'tone1', 'tone2', 'trace1', 'trace2'};
    % pre 
    phases{1}{1} = this_sess.csc.tvec(1);
    phases{1}{2} = this_sess.evts.t{2}(1);
    % post
    phases{2}{1} =  this_sess.evts.t{2}(2);
    phases{2}{2} =  this_sess.csc.tvec(end);
    %baseline
    tone_t_on = sort([this_sess.evts.t{4}(1:2:end); this_sess.evts.t{5}(1:2:end)]);
    phases{3}{1} = tone_t_on-30;
    phases{3}{2} = tone_t_on;

    % tone1 
    tone_t_1 = sort([this_sess.evts.t{4}(1:2:end)]) ;
    phases{4}{1} =  tone_t_1; 
    phases{4}{2} = tone_t_1 +20;

    % tone 2
    tone_t_2 = sort([this_sess.evts.t{5}(1:2:end)]);
    phases{5}{1} =  tone_t_2; 
    phases{5}{2} = tone_t_2 +20;

    % trace 1
    tone_t_1end = sort([this_sess.evts.t{4}(2:2:end)]);
    phases{6}{1} =  tone_t_1end;
    phases{6}{2} = tone_t_1end +20;

    % trace 2
    tone_t_2end = sort([this_sess.evts.t{5}(2:2:end)]);
    phases{7}{1} =  tone_t_2end;
    phases{7}{2} = tone_t_2end +20;

    for pp = 1:length(phases)
        % restrict data to the TOI
        this_ca1 = restrict(ca1_ts, phases{pp}{1}, phases{pp}{2});
        this_sub = restrict(sub_ts, phases{pp}{1}, phases{pp}{2});
        this_std = restrict(std_ts, phases{pp}{1}, phases{pp}{2});
        this_atyp = restrict(atyp_ts, phases{pp}{1}, phases{pp}{2});


        fprintf('Ca1 bursts in %s: <strong>%d</strong> \n',  phase_labels{pp}, length(this_ca1.t{1}))
        fprintf('Sub bursts in %s: <strong>%d</strong> \n',  phase_labels{pp}, length(this_sub.t{1}))
        fprintf('std Sub bursts in %s: <strong>%d</strong> \n',  phase_labels{pp}, length(this_std.t{1}))
        fprintf('atyp Sub bursts in %s: <strong>%d</strong> \n',  phase_labels{pp}, length(this_atyp.t{1}))

        % fill in the table. Coule be done in one line using a concatenation, but
        % this keeps it easy.
        mua_counts(pp,1) = length(this_ca1.t{1});
        mua_counts(pp,2) = length(this_sub.t{1});
        mua_counts(pp,3) = length(this_std.t{1});
        mua_counts(pp,4) = length(this_atyp.t{1});
        mua_phase_labels{pp} = phase_labels{pp}; 
    end
    % convert to proportions

    out.mua_prop{iS} = mua_counts ./ (mua_counts(:,3) + mua_counts(:,4));
    out.mua_rows{iS} = mua_phase_labels;

    clearvars -except iS out f_list  S_idx
end % end cross sessions.




%% intra session stats

% collect the CCF

mua_ccf = cell2mat(out.mua_ccf)';
mua_ccf_tvec = out.t_out{1};


% compile spike activity across atyp and std SWRs
Spk_std = cell2mat(out.Spk_sub_std);
Spk_atyp = cell2mat(out.Spk_sub_atyp);
Spk_loc = logical(cell2mat(out.Spk_loc));
Spk_pyr = cell2mat(out.Spk_rate) < 20;


% collect the behaviours.
hab_idx = contains(out.sess, 'TFCD1') | contains(out.sess, 'TFCD2');
train_idx = contains(out.sess, 'TFCD3');
test_idx = contains(out.sess, 'TFCD4') | contains(out.sess, 'TFCD5');
pox_idx = contains(out.sub, '3265'); 


% convert to a matrix
prop_mat = []; 
for ii = length(out.mua_prop):-1:1
    if pox_idx(ii)
    prop_mat(:,:,ii) = NaN(7,4); 

    else
    prop_mat(:,:,ii) = out.mua_prop{ii}; 
    end
end

mua_prop.all = prop_mat(:,:,~pox_idx )*100; 

mua_prop.hab = prop_mat(:,:,~pox_idx & hab_idx )*100; 
% mua_prop.train = prop_mat(:,:,~pox_idx & train_idx & k_idx); 
mua_prop.test = prop_mat(:,:,~pox_idx & test_idx )*100; 



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

%% plot the CCF for the Ca1 and Sub MUA peaks; 
    figure(1919)
    clf
    subplot(2,2,1)
    plot(mua_ccf_tvec*1000, mua_ccf./max(mua_ccf,[],2));
    hold on
    plot(mua_ccf_tvec*1000, mean(mua_ccf./max(mua_ccf,[],2)), 'k', 'LineWidth',2);
    
    legend('Sess 1', 'Sess 2', 'Sess 3', 'Sess 4', 'Sess 5', 'Mean', 'box', 'off')
    set(gca, 'xtick', -50:25:50, 'ytick', 0:.25:1)
    ylabel({'Cross-correlation'; '(normalized)'})
    xlabel('Time from Ca1 Burst (ms)')
    % axis square

    exportgraphics(gcf, [fig_dir filesep 'ccf_mua.pdf'], 'ContentType', 'vector');

%% get the proportion of SWRs per task phase

figure(3005)
clf

cond = 'all';
subplot(2,3,1)
cla
h1 = MS_bar_w_err(squeeze(mean(mua_prop.(cond)(1:2,3,:), 'omitmissing'))', squeeze(mean(mua_prop.(cond)(1:2,4,:), 'omitmissing'))', [c_ord(1,:); c_ord(1,:)], 1, 'ttest', 1:2, 1);
ylim([0 100])
set(gca, 'xtick', [1.5], 'XTickLabel', {'pre+post'})

legend(h1, 'Std', 'box', 'off')
title(['Sub SWRs in ' cond ' '])
ylabel('% of Bursts')

cond = 'all';
subplot(2,3,4)
cla
h1 = MS_bar_w_err(squeeze(mua_prop.(cond)(1,4,:))', squeeze(mua_prop.(cond)(2,4,:))', [c_ord(1,:)*.75; c_ord(1,:)], 1, 'ttest', 1:2);
h1.EdgeColor = 'none';
ylim([0 100])
set(gca, 'xtick', [1 2], 'XTickLabel', {'pre', 'post'})
ylabel('% of Atypical Bursts')
title(['Sub bursts in ' cond ' '])


% overall
cond = 'all';
subplot(2,3,2)
cla
MS_bar_w_err4(squeeze(mean(mua_prop.(cond)(1:2,4,:), 'omitmissing')), squeeze(mua_prop.(cond)(3,4,:)),...
    squeeze(mean(mua_prop.(cond)(4:5,4,:), 'omitmissing')), squeeze(mean(mua_prop.(cond)(6:7,4,:), 'omitmissing')),...
    c_ord(1:4,:), 1, 'anova2', 1:4, {'pre+post', 'baseline', 'tones', 'trace'});

ylim([0 100])
set(gca,'XTickLabel', {'pre+post', 'baseline', 'tones', 'trace'})

title(['Sub bursts in ' cond ' '])
ylabel('% of Atypical Bursts')

cond = 'hab';
subplot(2,3,3)
cla
MS_bar_w_err4(squeeze(mean(mua_prop.(cond)(1:2,4,:), 'omitmissing'))', squeeze(mua_prop.(cond)(3,4,:))',...
    squeeze(mean(mua_prop.(cond)(4:5,4,:), 'omitmissing'))', squeeze(mean(mua_prop.(cond)(6:7,4,:), 'omitmissing'))',...
    c_ord(1:4,:), 1, [], 1:4, {'pre+post', 'baseline', 'tones', 'trace'});

ylim([0 100])
set(gca,'XTickLabel', {'pre+post', 'baseline', 'tones', 'trace'})
title(['Sub bursts in ' cond ' '])


% cond = 'train';
% subplot(2,3,5)
% cla
% MS_bar_w_err4(squeeze(mean(mua_prop.(cond)(1:2,3,:), 'omitmissing'))', squeeze(mean(mua_prop.(cond)(3,3,:), 'omitmissing'))',...
%     squeeze(mean(mua_prop.(cond)(4:5,3,:), 'omitmissing'))', squeeze(mean(mua_prop.(cond)(6:7,3,:), 'omitmissing'))',...
%     c_ord(1:4,:), 1, [], 1:4, {'pre+post', 'baseline', 'tones', 'trace'})
% 
% ylim([0 1.25])
% set(gca,'XTickLabel', {'pre+post', 'baseline', 'tones', 'trace'})
% title(['Sub SWRs in ' cond ' '])

cond = 'test';
subplot(2,3,6)
cla
MS_bar_w_err4(squeeze(mean(mua_prop.(cond)(1:2,4,:), 'omitmissing')),squeeze(mua_prop.(cond)(3,4,:)),...
    squeeze(mean(mua_prop.(cond)(4:5,4,:), 'omitmissing')), squeeze(mean(mua_prop.(cond)(6:7,4,:), 'omitmissing')),...
    c_ord(1:4,:), 1, 'anova2', 1:4, {'pre+post', 'baseline', 'tones', 'trace'});

ylim([0 100])
set(gca,'XTickLabel', {'pre+post', 'baseline', 'tones', 'trace'})

title(['Sub bursts in ' cond ' '])

    exportgraphics(gcf, [fig_dir filesep 'task_phase_mua.pdf'], 'ContentType', 'vector');



