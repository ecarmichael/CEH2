%% sandbox_NAMI_Sub_SWR
% 

%pox2217_TFCD1
% csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Pox/Pox2217_2026-06-16_16-29-44_TFC_D1/Record Node 117';
% csc_idx = 1:4:96;
% ts_prime = 0;
% csc_idx = {'CH5', 'CH9','CH63' 'CH51' 'CH55', 'CH63'}; % no good Sub SWR


%pox2217_TFCD1
csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Pox/Pox2217_2026-06-10_11-27-27_LT2/Record Node 117';
csc_idx = 1:4:96;
ts_prime = 0;
% csc_idx = {'CH5', 'CH9','CH63' 'CH51' 'CH55', 'CH63'}; % no good Sub SWR

%pox3567_TFCD1
% csc_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\Pox\Pox3567_2026-06-20_17-44-01_TFCD1\Record Node 117';
% % csc_idx = 1:4:96;
% ts_prime = 0;
% csc_idx = {'CH51', 'CH143'};


%pox3256_TFCD1  %% no ripples at all? 
% csc_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\Pox\Pox3265_2026-06-16_18-09-08_TFC_D1\Record Node 117';
% csc_idx = 1:4:96;
% ts_prime = 0;
% csc_idx = {'CH51', 'CH143'};

%pox3256_TFCD3  %% no ripples at all? 
% csc_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\Pox\Pox3265_2026-06-16_18-09-08_TFC_D1\Record Node 117';
% csc_idx = 1:4:96;
% ts_prime = 0;
% csc_idx = {'CH51', 'CH143'};
%% load the spikes if present
%
% params = OE_load_params(phy_dir);
%
% data.S = OE_phy2TS(phy_dir, params);
%
% ts_prime = readNPY([phy_dir filesep 'timestamps.npy']);
%% load the csc


if isempty(csc_dir)
    csc = [];
    OE_evts =[];
else
    csc_list = dir([csc_dir filesep '*CH*.continuous']);

    % sort the csc based on channel number.
    for ii = length(csc_list):-1:1
        ch_idx =  strfind(csc_list(ii).name,'_CH');
        con_idx =  strfind(csc_list(ii).name,'.continuous');
        csc_num(ii) = str2double(csc_list(ii).name(ch_idx+3:con_idx));
        csc_name{ii} = csc_list(ii).name;
    end

    % sort
    [~, sort_idx] = sort(csc_num);
    csc_list = csc_list(sort_idx);

    if ~isempty(csc_idx)  && isnumeric(csc_idx(1))
        csc_list(~ismember(1:length(csc_list), csc_idx))= [];

    elseif  ~isempty(csc_idx)  && iscell(csc_idx) % if csc_idx is a string look for those patterns

        temp_list = [];

        for ii = 1:length(csc_idx)
            temp_idx = find(contains(csc_name, [csc_idx{ii} '.continuous']));
            if length(temp_idx) > 1
                error('Too many csc files containing this name')
            end
            if ~isempty(temp_idx)
                temp_list(ii).name= csc_name{temp_idx};
                temp_list(ii).folder= csc_list(1).folder;
            end
        end

        csc_list = temp_list;

    end


    csc= []; labels = [];

    % get the first channel with the time vector.
    ii = 1;
    [data, tvec, info] = load_open_ephys_data([csc_list(ii).folder filesep csc_list(ii).name]);
    csc = tsd(tvec, data);
    labels{ii} = info.header.channel;
    csc.cfg.hdr{ii} = info.header;
    csc.cfg.hdr{ii}.SamplingFrequency = info.header.sampleRate;
    % csc.data = [csc.data, NaN(length(csc.data),length(csc_list)-1)]; % pad the NaNa

    for ii = length(csc_list):-1:2
        [data, ~, info] = load_open_ephys_data([csc_list(ii).folder filesep csc_list(ii).name]);
        csc.data(:,ii) =  data;
        labels{ii} = info.header.channel;
        csc.cfg.hdr{ii} = info.header;
        csc.cfg.hdr{ii}.SamplingFrequency = info.header.sampleRate;

    end
    fs = csc.cfg.hdr{1}.SamplingFrequency;

    csc.data = csc.data';
    csc.label = labels;


    % csc.tvec = ;
    % csc.tvec = csc.tvec - csc.tvec(1) + (csc.tvec(1) - ts_prime(1)); % zero out the csc.
    csc.tvec = csc.tvec - ts_prime(1); % zero out the csc.

    cfg_in.decimateFactor = 15;
    csc = decimate_tsd(cfg_in, csc);

    % load the OE version of the events.
    evts_list = dir([csc_dir filesep '*Data*.events']);

    OE_evts = OE_LoadEvents([evts_list.folder filesep evts_list.name], fs);

end



%% csc check

figure(1010)
clf;
hold on
for ii = 1:size(csc.data,1)

    plot(csc.tvec, csc.data(ii,:)+ii*500);
    lab{ii} = csc.label{ii};
    y_t(ii) = median(csc.data(ii,:)+ii*500);
end
set(gca, 'YTick', y_t, 'YTickLabel', lab)

% vline(SWR_evts.t{2}, 'r')

% %% grab the SWR times
% evts_list = dir([swr_dir filesep '*Data*.events']);
% SWR_evts = OE_LoadEvents([evts_list.folder filesep evts_list.name], fs);
%
% iRi = diff(SWR_evts.t{2});
% keep_idx = iRi <.05;
%
% SWR_evts.t{2}(keep_idx) = [];

%% get teh movement if present in the evts

move_ts = ts({OE_evts.t{contains(OE_evts.label, '8')}});

mov_rate  = MS_spike2rate(move_ts, csc.tvec);


cfg_mov = [];
cfg_mov.threshold = .001;
cfg_mov.dcn = '<';
cfg_mov.operation = '<';
cfg_mov.minlen = .5;

mov_iv = TSDtoIV(cfg_mov, mov_rate);

% pad movement

cfg_resize.d = [+.5 -.5];

mov_iv = ResizeIV(cfg_resize, mov_iv);


csc_r = restrict(csc, mov_iv);

%%   Detect CA1 SWRS

    swrs_ca1 = MS_SWR_detector(csc_r,csc.label{1});


%%   Detect Sub SWRS
close all
    swrs_sub = MS_SWR_detector(csc_r,csc.label{2});


%% collect the data
load("all_TFC.mat")
this_name = 'pox_3567_TFCD1';

all_TFC.(this_name).csc = csc;

all_TFC.(this_name).swrs_ca1= swrs_ca1;
all_TFC.(this_name).swrs_sub= swrs_sub;


all_TFC.(this_name).mov_iv= mov_iv;
all_TFC.(this_name).mov_vec= mov_rate;

OE_evts.label = {'context' 'TFC_on' 'mov', 'tone1', 'tone2', 'lick', 'eye'}; 

all_TFC.(this_name).evts=OE_evts;

save('all_TFC.mat', 'all_TFC')



%% chopping data

% grab some data
this_data = all_TFC.pox_3567_TFCD1;

% restrict data to the pre task baseline
% ca1
swr_ca1_pre = restrict(this_data.swrs_ca1, this_data.csc.tvec(1), this_data.evts.t{2}(1)); 

%sub 
swr_sub_pre = restrict(this_data.swrs_sub, this_data.csc.tvec(1), this_data.evts.t{2}(1)); 

fprintf('Ca1 SWRs in Pre: <strong>%d</strong> \n', length(swr_ca1_pre.tstart))

fprintf('Sub SWRs in Pre: <strong>%d</strong> \n', length(swr_sub_pre.tstart))


% restrict data to the post period. 

% ca1
swr_ca1_post = restrict(this_data.swrs_ca1, this_data.evts.t{2}(2), this_data.csc.tvec(end)); 

%sub 
swr_sub_post = restrict(this_data.swrs_sub, this_data.evts.t{2}(2), this_data.csc.tvec(end)); 

fprintf('Ca1 SWRs in Post: <strong>%d</strong> \n', length(swr_ca1_post.tstart))

fprintf('Sub SWRs in Post: <strong>%d</strong> \n', length(swr_sub_post.tstart))


% get the baseline periods (30s before any tone) inside the TFC
dur = 30; 
tone_t_on = sort([this_data.evts.t{4}(1:2:end); this_data.evts.t{5}(1:2:end)]);

swr_ca1_base = restrict(this_data.swrs_ca1, tone_t_on-30, tone_t_on); 
swr_sub_base = restrict(this_data.swrs_sub, tone_t_on-30, tone_t_on); 

fprintf('Ca1 SWRs in baseline: <strong>%d</strong>   |    Sub SWRs in baseline: <strong>%d</strong> \n',...
    length(swr_ca1_base.tstart),  length(swr_sub_base.tstart))


% get the tone 1 periods 




% get the tone 2 periods



% tone 1 trace period (20 sec after end of tone)



% tone 2 trace period 








%% plots to check everything

figure(101)  % make a plot (with specific number)
clf % clear the current plot

hold on % keep things on the same plot

% plot the data
plot(all_TFC.pox_3567_TFCD1.csc.tvec, all_TFC.pox_3567_TFCD1.csc.data(1:2,:))


% add some vertical lines
xline(tone_t_on-30, 'y')

% add some fancy verical lines
xline(all_TFC.pox_3567_TFCD1.evts.t{5}(1:2:end), 'k', 'LineWidth',3)










%% LATER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(998)
clf
set(gcf,'Units','inch','position',[3 3 3.5 2.5]);
[~, ~, ~, p, stats] =MS_bar_w_err(all_ctr_type(:,1), all_pox_type(:,1), c_ord, 1, 'ttest2', [1 1.5]);

% [~, ~, ~, p, stats] = MS_bar_w_err(all_ctr_dur, all_pox_dur, c_ord, 1, 'ttest2', [1 1.5]);
% sc{1}.MarkerFaceColor = c_ord(1,:); sc{2}.MarkerFaceColor = c_ord(2,:);
fprintf('SWR rate: t(%.0f): %.2f, p = %.3f\n',stats.df, stats.tstat, p )

xlim([.75 4.75]); ylim([0 .3])
set(gca, "XTick", [1.25 2.25 3.25 4.25], 'XTickLabel', {'1x', '2x', '3x', 'Other'},...
    'YTick', 0:.1:.3, 'TickDir', 'out')
ylabel('SWR rate (Hz)')
vline(5, 'k')


figure(997)
clf
set(gcf,'Units','inch','position',[3 3 3.5 2.5]);
% [~, ~, ~, p, stats] = MS_bar_w_err(all_ctr_dur*1000, all_pox_dur*1000, c_ord, 1, 'ttest2', [1 1.5]);
MS_rain_plot([all_ctr_dur; all_pox_dur]', [zeros(1,length(all_ctr_dur)), ones(1,length(all_pox_dur))]', c_ord, 'ttest2', [1 3]);
% [~, ~, ~, p, stats] = MS_bar_w_err(all_ctr_dur*1000, all_pox_dur*1000, c_ord, 1, 'ttest2', [1 1.5]);
% sc{1}.MarkerFaceColor = c_ord(1,:); sc{2}.MarkerFaceColor = c_ord(2,:);
fprintf('SWR dur: t(%.0f): %.2f, p = %.3f\n',stats.df, stats.tstat, p )

ylim([.75 6.75]); xlim([0 .12])
set(gca, "yTick", [1.25 3.25 ], 'yTickLabel', {'aTau-', 'aTau+'}, 'TickDir', 'out', 'XTick', 0:.05:.1);
set(gca, 'XTickLabel', get(gca, 'xtick').*1000)
%'YTick', 0:.1:.3)
xlabel('SWR duration (ms)')
%% quick glm

t_type = table([1:(size(all_pox_prob,1)+size(all_ctr_prob,1))]', [repmat({'Pox'}, size(all_pox_type,1), 1); repmat({'Ctrl'}, size(all_ctr_type,1), 1)], [all_pox_type(:,1); all_ctr_type(:,1)],...
    [all_pox_type(:,2); all_ctr_type(:,2)], [all_pox_type(:,2); all_ctr_type(:,2)],[all_pox_type(:,4); all_ctr_type(:,4)],...
    'VariableNames',{'Subject', 'Group', 'SWR1', 'SWR2', 'SWR3', 'SWR4'});

withinDesign = table([1 2 3 4]', 'VariableNames', {'type'});

rm = fitrm(t_type, 'SWR1-SWR4~ Group', 'WithinDesign', withinDesign);

ranovaTable = ranova(rm);

fprintf('SWR type\n')
display(ranovaTable)

Tukey_table = multcompare(rm,'Group', 'By', 'type');




t_prob = table([1:(size(all_pox_prob,1)+size(all_ctr_prob,1))]', [repmat({'Pox'}, size(all_pox_prob,1), 1); repmat({'Ctrl'}, size(all_ctr_prob,1), 1)], [all_pox_prob(:,1); all_ctr_prob(:,1)],...
    [all_pox_prob(:,2); all_ctr_prob(:,2)], [all_pox_prob(:,2); all_ctr_prob(:,2)],[all_pox_prob(:,4); all_ctr_prob(:,4)],...
    'VariableNames',{'Subject', 'Group', 'SWR1', 'SWR2', 'SWR3', 'SWR4'});

withinDesign = table([1 2 3 4]', 'VariableNames', {'type'});

rm = fitrm(t_prob, 'SWR1-SWR4~ Group', 'WithinDesign', withinDesign);

ranovaTable = ranova(rm);

fprintf('SWR Prob\n')
display(ranovaTable)

Tukey_table = multcompare(rm,'Group', 'By', 'type');


%% simple center offset histogram

sub_cent = IVcenters(sub_swrs);
ca1_cent = IVcenters(ca1_swrs);
pair_swr = [];
for ii = length(ca1_cent):-1:1
    this_d = sub_cent - ca1_cent(ii);


    co_evt = find(abs(this_d) < .50);

    if isempty(co_evt) || length(co_evt) > 1
        pair_swr(ii) = NaN;
    else
        pair_swr(ii) = this_d(co_evt);
    end
end

figure(8898)
clf

histogram(pair_swr(~isnan(pair_swr))*1000, [-25:1:25])
xline(0, '--k', 'LineWidth',1)
axis('square')