%% sandbox HF_SWR_spec
% pox1_2
% csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_14-03-46_SWR/Record Node 117'; 
% swr_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_14-03-46_SWR/Record Node 143'; 
% csc_idx = [24, 34];

% % pox 1_2
% phy_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_14-03-46_SWR/kilosort4'; 
% csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_14-03-46_SWR/Record Node 117'; 
% swr_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_14-03-46_SWR/Record Node 143'; 
% csc_idx = [ 15    26     28    30    32    34]; % 26

% pox 2_1
% phy_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-16-12_SWR/kilosort4'; 
% csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-16-12_SWR/Record Node 117'; 
% swr_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-16-12_SWR/Record Node 143'; 
% % csc_idx = [ 15    26     28    30    32    34]; % 26

% % % % pox_2_2
% % % phy_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-30-53_SWR/kilosort4'; 
% % % csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-30-53_SWR/Record Node 117'; 
% % % swr_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-30-53_SWR/Record Node 143'; 
% % % % csc_idx = 26; %[ 5:16    20:30    32    33:36]; % 26

% % % % pox_2_3
% % % phy_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-39-33_SWR2/kilosort4'; 
% % % csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-39-33_SWR2/Record Node 117'; 
% % % swr_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-39-33_SWR2/Record Node 143'; 
% % % csc_idx = [24 28 30]; %26; %[ 5:16    20:30    32    33:36]; % 26

% % pox_2_4
% phy_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-53-28_SWR2/kilosort4'; 
% csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-53-28_SWR2/Record Node 117'; 
% swr_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-53-28_SWR2/Record Node 143'; 
% csc_idx = 26; %[ 5:16    20:30    32    33:36]; % 26

% ctrl
% csc_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\test_data\H1b32026-03-14_22-43-13_dSub_SWR\Record Node 117'

% ctrl H1
% csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/H1b32026-03-14_22-43-13_dSub_SWR/Record Node 117'
% csc_idx = [24 28 30]; 

% ctrl 3b3
% csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/TFC/H3b3/HF3b3_2026-03-03_13-21-54_D2_Hab/Record Node 117'
% csc_idx = [14]; 

% ctrl_h1b3
% Xcsc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/TFC/H1b3/habituation/HF1b3_2026-03-01_18-53-32_D1_Hab/Record Node 117'; 
% csc_idx = [29:2:42]; 

% ctrl_h2b3_hab 2
csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/TFC/H2b3/HF2b3_2026-03-02_13-26-03_D1_Hab/Record Node 143'; 
% csc_idx = [ 1:2:34]; 
ts_prime = 0; 

%% load the spikes if present

params = OE_load_params(phy_dir);

data.S = OE_phy2TS(phy_dir, params);

ts_prime = readNPY([phy_dir filesep 'timestamps.npy']);
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
for ii = 1:length(csc_list)

    if ii == 1
        [data, tvec, info] = load_open_ephys_data([csc_list(ii).folder filesep csc_list(ii).name]);
        csc = tsd(tvec, data);
        labels{ii} = info.header.channel;
        csc.cfg.hdr{ii} = info.header;
        csc.cfg.hdr{ii}.SamplingFrequency = info.header.sampleRate;

    else
        [data, ~, info] = load_open_ephys_data([csc_list(ii).folder filesep csc_list(ii).name]);
        csc.data =[csc.data, data];
        labels{ii} = info.header.channel;
        csc.cfg.hdr{ii} = info.header;
        csc.cfg.hdr{ii}.SamplingFrequency = info.header.sampleRate;
    end
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

%


if length(csc_idx) > 2
csc.data(end+1,:) = mean(csc.data); 

csc.label{end+1} = 'mean'; 
csc.cfg.hdr{end+1} = csc.cfg.hdr{end}; 
end


%% grab the SWR times
evts_list = dir([swr_dir filesep '*Data*.events']);
SWR_evts = OE_LoadEvents([evts_list.folder filesep evts_list.name], fs);

iRi = diff(SWR_evts.t{2}); 
keep_idx = iRi <.05; 

SWR_evts.t{2}(keep_idx) = []; 
%% trim if needed
% ctrl_h1
% swr_cut_iv = iv([csc.tvec(1) 385 460 940 1580], [350  440 560 1140 1840 ])

%ctrl 3b3
swr_cut_iv = iv([360 1340 1520 3180], [500 1440 1580 3240]); 

% pox1_1
% swr_cur_iv = iv([csc.tvec(1) 555], [515 csc.tvec(end)]); 
% theta_cur_iv = iv([3400 3706],  [3460  3805]); 

% % pox2_4
% swr_cur_iv = iv([csc.tvec(1) 3460 3484.5 3513.35 3805], [3400 3483.5 3512.35 3706 csc.tvec(end)]); 
% theta_cur_iv = iv([3400 3706],  [3460  3805]); 

% csc_theta = restrict(csc, theta_cur_iv); 
%%

swrs = MS_SWR_detector(csc, 'CH8');

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

%% get Single Double Tripples as per Yamamato & Tonegawa Neuron 2017
% https://www.cell.com/neuron/fulltext/S0896-6273(17)30857-7#sec-4
% this_swr = ctrl_h3b3.swrs; 
% this_csc = ctrl_h3b3.csc;  thi

swr_type  = NaN(size(this_swr.tstart));

for ii = 1:length(this_swr.tstart)-1
    d = this_swr.tstart(ii+1:end) - this_swr.tend(ii); 

    if (length(d) > 2) &&(d(1) <= .200) && (d(1) >.070) && (d(2) <= .400) && (d(2) >.070)
    swr_type(ii) = 3;
    elseif (d(1) <= .200) && (d(1) >.070)
        swr_type(ii) = 2;
    elseif d(1) > .200 && (d(1) >.070)
        swr_type(ii) = 1;
    else
    swr_type(ii) = 4;
    end
end
swr_type(end) = []; 

t_idx = [swr_type(1:end-2) == 3; 0; 0]; 
swr_type([find(t_idx)+1 find(t_idx)+2]) = []; 

d_idx = [swr_type(1:end-1) == 2; 0]; 
swr_type(find(d_idx)+1) = []; 

figure(101)
clf
histogram(swr_type, 'Normalization','probability')
set(gcf,'Units','inch','position',[3 3 5 4]);
ylabel('probability'); xlabel('SWR type')

swr_types_rate = [sum(swr_type == 1) ./ (length(this_csc.tvec) / this_csc.cfg.hdr{1}.SamplingFrequency),...
    sum(swr_type == 2) ./ (length(this_csc.tvec) / this_csc.cfg.hdr{1}.SamplingFrequency), ...
        sum(swr_type == 3) ./ (length(this_csc.tvec) / this_csc.cfg.hdr{1}.SamplingFrequency),...
            sum(swr_type == 4) ./ (length(this_csc.tvec) / this_csc.cfg.hdr{1}.SamplingFrequency)];

swr_s = SelectIV([], this_swr, swr_type == 1); 
swr_d = SelectIV([], this_swr, swr_type == 2); 
swr_t = SelectIV([], this_swr, swr_type == 3); 

%% collect the data
load("all_data.mat")
this_name = 'ctrl_h2b3_d3'; 

all_data.(this_name).csc = csc; 

all_data.(this_name).swrs= swrs;
all_data.(this_name).swr_rate = length(swrs.tstart) / (length(csc.tvec) ./ csc.cfg.hdr{1}.SamplingFrequency); 
all_data.(this_name).swr_type= swr_type;
all_data.(this_name).swr_type_rate= swr_types_rate;

all_data.(this_name).swr_s= swr_s; 
all_data.(this_name).swr_d= swr_d; 
all_data.(this_name).swr_t= swr_t; 

save('all_data.mat', 'all_data')
%% histo of SWR diffs

figure(1919)
subplot(2,2,1)
histogram(diff(SWR_evts.t{2}), 0:.01:1)
vline([110, 220, 330]./1000)

subplot(2,2,3)
% pie([])

subplot(2,2,2)
histogram(diff(IVcenters(swrs)), 0:.01:1)
vline([110, 220, 330]./1000)



%% get the spectrogram


%%

this_csc = ctrl_h1.csc;
this_swr = ctrl_h1.swrs; 

csc_idx = contains(this_csc.label, 'mean');

csc_ft = this_csc; 
csc_ft.data = this_csc.data(find(csc_idx),:); 
csc_ft.label = []; 
csc_ft.label{1} = this_csc.label{find(csc_idx)}; 



swr_centers = IVcenters(this_swr); 

[csc_ft_out, TFR] = Triggered_Spec_FT(csc_ft, swr_centers, [], 80:.5:200, [], [-.5 .5]);
z_idx = nearest_idx3(0, TFR.time); 
b_idx = [nearest_idx3(-.5, TFR.time), nearest_idx3(-.2, TFR.time)];  

set(gcf,'Units','inch','position',[3 3 5 4]);
xlim([-.05 .05]); ylabel('frequency (hz)'); xlabel('time from ripple center (ms)')
set(gca,'xtick', -.05:0.025:.05); 
set(gca, 'XTickLabel', get(gca, 'XTick')*1000)

for ii = length(swr_centers):-1:1
    [~, TFR] = Triggered_Spec_FT(csc_ft, swr_centers(ii), [], 80:.5:200,  [], [-.5 .5], 0);
    psds(ii, :) = (squeeze(TFR.powspctrm(1,:,z_idx)));% - mean(squeeze(TFR.powspctrm(1,:,b_idx(1):b_idx(2))), 2, 'omitmissing')')./std(mean(squeeze(TFR.powspctrm(1,:,b_idx(1):b_idx(2))), 2, 'omitmissing')');
end


%% average PSD per ripple

% this_psd = squeeze(TFR.powspctrm(1,:,z_idx)); 

figure(1013)
clf
hold on

yMean = mean(ctrl_h2.psds, 1);
ySEM = std(ctrl_h2.psds, 0, 1) / sqrt(size(ctrl_h2.psds, 1));

shadedErrorBar(ctrl_h2.TFR.freq, yMean, ySEM, 'lineprops', '-k');


yMean = mean(pox1.psds, 1);
ySEM = std(pox1.psds, 0, 1) / sqrt(size(pox1.psds, 1));

shadedErrorBar(pox1.TFR.freq, yMean, ySEM, 'lineprops', '-r');

%% theta gamma mod

[CoMo,Phi_f, Amp_f]= MS_phase_freq([], csc_theta, [4, 12], [20 120]);

figure(10111)
clf
imagesc(Phi_f, Amp_f, CoMo')
title('Pox2_4')
%% spectrogram 


%%  bar plot with some stats

all_pox_rate = [all_data.pox1_1.swr_rate, all_data.pox1_2.swr_rate , all_data.pox2_4.swr_rate]; 
all_ctr_rate = [all_data.ctrl_h1.swr_rate, all_data.ctrl_h2.swr_rate , all_data.ctrl_h3b3.swr_rate, all_data.ctrl_h2b3_d1.swr_rate,...
    all_data.ctrl_h2b3_d3.swr_rate]; 

% type rate
all_pox_type = [all_data.pox1_1.swr_type_rate;   all_data.pox1_2.swr_type_rate; all_data.pox2_4.swr_type_rate]; 
all_ctr_type = [all_data.ctrl_h1.swr_type_rate; all_data.ctrl_h2.swr_type_rate; all_data.ctrl_h3b3.swr_type_rate;...
    all_data.ctrl_h2b3_d1.swr_type_rate; all_data.ctrl_h2b3_d3.swr_type_rate]; 


% get the ripple durations
f = fieldnames(all_data); 
for ii = 1:length(f)
    [all_data.(f{ii}).swr_type_prob , b] = histcounts(all_data.(f{ii}).swr_type,.5:1:4.5, 'Normalization','probability'); 

    % type rate
    all_data.(f{ii}).swr_dur = all_data.(f{ii}).swrs.tend - all_data.(f{ii}).swrs.tstart;

end

% type prob
all_pox_prob = [all_data.pox1_1.swr_type_prob;   all_data.pox1_2.swr_type_prob; all_data.pox2_4.swr_type_prob]; 
all_ctr_prob = [all_data.ctrl_h1.swr_type_prob; all_data.ctrl_h2.swr_type_prob; all_data.ctrl_h3b3.swr_type_prob;...
    all_data.ctrl_h2b3_d1.swr_type_prob; all_data.ctrl_h2b3_d3.swr_type_prob]; 



all_pox_dur = [all_data.pox1_1.swr_dur; all_data.pox1_2.swr_dur; all_data.pox2_4.swr_dur]; 
all_ctr_dur = [all_data.ctrl_h1.swr_dur; all_data.ctrl_h2.swr_dur; all_data.ctrl_h3b3.swr_dur; all_data.ctrl_h2b3_d1.swr_dur; ...
    all_data.ctrl_h2b3_d3.swr_dur]; 

% plots
c_ord = [0.502 0.502 0.502 ; .8 0 0]; 


figure(999)
set(gcf,'Units','inch','position',[3 3 3.5 2.5]);
clf; 
hold on

[h, eb, sc] = MS_bar_w_err(all_ctr_prob(:,1), all_pox_prob(:,1), c_ord, 1, 'ttest2', [1 1.5])
% sc{1}.MarkerFaceColor = c_ord(1,:); sc{2}.MarkerFaceColor = c_ord(2,:);

MS_bar_w_err(all_ctr_prob(:,2), all_pox_prob(:,2), c_ord, 1, 'ttest2', [2 2.5])
MS_bar_w_err(all_ctr_prob(:,3), all_pox_prob(:,3), c_ord, 1, 'ttest2', [3 3.5])
MS_bar_w_err(all_ctr_prob(:,4), all_pox_prob(:,4), c_ord, 1, 'ttest2', [4 4.5])

xlim([.75 4.75]); %ylim([0 .3])
set(gca, "XTick", [1.25 2.25 3.25 4.25], 'XTickLabel', {'1x', '2x', '3x', 'Other'},...
    'YTick', 0:.25:1, 'TickDir', 'out')
ylabel('probability')
vline(5, 'k')
% legend('aTau-', 'aTau+')


figure(909)
set(gcf,'Units','inch','position',[3 3 3.5 2.5]);
clf
hold on

MS_bar_w_err(all_ctr_type(:,1), all_pox_type(:,1), c_ord, 1, 'ttest2', [1 1.5]); 
% sc{1}.MarkerFaceColor = c_ord(1,:); sc{2}.MarkerFaceColor = c_ord(2,:);

MS_bar_w_err(all_ctr_type(:,2), all_pox_type(:,2), c_ord, 1, 'ttest2', [2 2.5])
MS_bar_w_err(all_ctr_type(:,3), all_pox_type(:,3), c_ord, 1, 'ttest2', [3 3.5])
MS_bar_w_err(all_ctr_type(:,4), all_pox_type(:,4), c_ord, 1, 'ttest2', [4 4.5])

xlim([.75 4.75]); ylim([0 .3])
set(gca, "XTick", [1.25 2.25 3.25 4.25], 'XTickLabel', {'1x', '2x', '3x', 'Other'},...
    'YTick', 0:.1:.3, 'TickDir', 'out');
ylabel('SWR rate (Hz)')
vline(5, 'k')
% legend('aTau-', 'aTau+')



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
MS_rain_plot([all_ctr_dur; all_pox_dur]', [zeros(1,length(all_ctr_dur)), ones(1,length(all_pox_dur))]', c_ord, 'ttest2', [1 3])
% [~, ~, ~, p, stats] = MS_bar_w_err(all_ctr_dur*1000, all_pox_dur*1000, c_ord, 1, 'ttest2', [1 1.5]);
% sc{1}.MarkerFaceColor = c_ord(1,:); sc{2}.MarkerFaceColor = c_ord(2,:);
fprintf('SWR dur: t(%.0f): %.2f, p = %.3f\n',stats.df, stats.tstat, p )

ylim([.75 6.75]); xlim([0 .12])
set(gca, "yTick", [1.25 3.25 ], 'yTickLabel', {'aTau-', 'aTau+'}, 'TickDir', 'out', 'XTick', [0:.05:.1]);
set(gca, 'XTickLabel', get(gca, 'xtick').*1000)
    %'YTick', 0:.1:.3)
xlabel('SWR duration (ms)')
%% quick glm

t_type = table([1 2 3 4 5 6 7 8]', [repmat({'Pox'}, 3, 1); repmat({'Ctrl'}, 5, 1)], [all_pox_type(:,1); all_ctr_type(:,1)],...
    [all_pox_type(:,2); all_ctr_type(:,2)], [all_pox_type(:,2); all_ctr_type(:,2)],[all_pox_type(:,4); all_ctr_type(:,4)],...
    'VariableNames',{'Subject', 'Group', 'SWR1', 'SWR2', 'SWR3', 'SWR4'}); 

withinDesign = table([1 2 3 4]', 'VariableNames', {'type'});

rm = fitrm(t_type, 'SWR1-SWR4~ Group', 'WithinDesign', withinDesign);

ranovaTable = ranova(rm); 

fprintf('SWR type\n')
display(ranovaTable)

Tukey_table = multcompare(rm,'Group', 'By', 'type'); 




t_prob = table([1 2 3 4 5 6 7 8]', [repmat({'Pox'}, 3, 1); repmat({'Ctrl'}, 5, 1)], [all_pox_prob(:,1); all_ctr_prob(:,1)],...
    [all_pox_prob(:,2); all_ctr_prob(:,2)], [all_pox_prob(:,2); all_ctr_prob(:,2)],[all_pox_prob(:,4); all_ctr_prob(:,4)],...
    'VariableNames',{'Subject', 'Group', 'SWR1', 'SWR2', 'SWR3', 'SWR4'}); 

withinDesign = table([1 2 3 4]', 'VariableNames', {'type'});

rm = fitrm(t_prob, 'SWR1-SWR4~ Group', 'WithinDesign', withinDesign);

ranovaTable = ranova(rm); 

fprintf('SWR Prob\n')
display(ranovaTable)

Tukey_table = multcompare(rm,'Group', 'By', 'type'); 


