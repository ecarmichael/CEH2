%% sandbox_HF_OE_Peth

%% HF2b2_D1_red_resp [only 2 cells?]
evts_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D1/HF2b2_2025-12-16_11-18-44_Red_respose_Ca1/Record Node 112/experiment1/recording1/events/Intan_RHD_USB-108.Rhythm Data/TTL' ;
csc_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D1/HF2b2_2025-12-16_11-18-44_Red_respose_Ca1/Record Node 117'; 
phy_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Kilo_inter/HF2b2_red_response';
vr_fname = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D1/HF2b2_D1_2025-12-16_12-47-39.csv'; 
csc_idx = [3 11]; 
swr_ch = 1; 
save_name = 'HF2b2_D1_red_resp'; 
TTL = {'6', '7'};

%% HF2b2_D1
evts_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D1/HF2b2_2025-12-16_12-39-09_D1/Record Node 112/experiment2/recording1/events/Intan_RHD_USB-108.Rhythm Data/TTL' ;
csc_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D1/HF2b2_2025-12-16_12-39-09_D1/Record Node 117'; 
phy_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Kilo_inter/HF2b2_D1';
vr_fname = []; 
csc_idx = [3 11]; 
swr_ch = 1; 
save_name = 'HF2b2_D1'; 
TTL = {'6', '7'};
    %% HF2b2_D3
evts_dir = ('/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D5/HF2b2_2025-12-18_13-47-24_D3_2_opto_only/Record Node 112/experiment1/recording1/events/Intan_RHD_USB-108.Rhythm Data/TTL') ;
csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D3/HF2b2_2025-12-18_13-47-24_D3_2_opto_only/Record Node 117';
phy_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Kilo_inter/HF2b2_D3';
vr_fname = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D3/HF2b2_D3_2025-12-18_14-08-34.csv';


%% HF2b2_D5
evts_dir = ('/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D5/HF2b2_2026-01-02_13-50-31_D5/Record Node 112/experiment1/recording1/events/Intan_RHD_USB-108.Rhythm Data/TTL') ;
csc_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D5/HF2b2_2026-01-02_13-50-31_D5/Record Node 117';
phy_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Kilo_inter/HF2b2_D5';
vr_fname = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D5/h2b2_VR_D5_2026-01-02_14-09-33.csv';
csc_idx = [1 6 11];
swr_ch = 2;
TTL = {'6', '7'};
save_name = 'HF2b2_D5';
%% HF2b2_D5_opto_1
% evts_dir = ('/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/HF2b2_2026-01-02_13-15-51_SS_test/Record Node 112/experiment1/recording1/events/Intan_RHD_USB-108.Rhythm Data/TTL') ;
% csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/HF2b2_2026-01-02_13-15-51_SS_test/Record Node 117';
% phy_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Kilo_inter/HF2b2_D5_Opto1';
% vr_fname = []; %'/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D3/HF2b2_D3_2025-12-18_14-08-34.csv';
% csc_idx = [1, 4, 6];
% swr_ch = 2;
%
% save_name = 'HF2b2_D5_opto1';
%% HF2b2_D5_opto_2
evts_dir = ('/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/HF2b2_2026-01-02_13-21-35_SS_test2/Record Node 112/experiment1/recording1/events/Intan_RHD_USB-108.Rhythm Data/TTL') ;
csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/HF2b2_2026-01-02_13-21-35_SS_test2/Record Node 117';
phy_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Kilo_inter/HF2b2_SS_D5_Opto_cells';
vr_fname = []; %'/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D3/HF2b2_D3_2025-12-18_14-08-34.csv';
csc_idx = [1, 4, 6];
swr_ch = 2;
save_name = 'HF2b2_D5_opto2';



%% HF2b2_D6_opto_1
evts_dir = '\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\GoNoGo\HF2b2_D6\HF2b2_2026-01-03_21-22-41_D6_opto_1\Record Node 112\experiment1\recording1\events\Intan_RHD_USB-108.Rhythm Data\TTL' ;
csc_dir = '\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\GoNoGo\HF2b2_D6\HF2b2_2026-01-03_21-22-41_D6_opto_1\Record Node 117'; 
phy_dir = '\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\Kilo_inter\HF2b2_D6_opto';
vr_fname = []; %'/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D3/HF2b2_D3_2025-12-18_14-08-34.csv'; 
csc_idx = [3 7 8]; 
swr_ch = 2; 
save_name = 'HF2b2_D6_opto1'; 

%% SOM2_Opto
evts_dir = '\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\test_data\SOM2_2025-10-22_15-21-16_opto_test_wide_probe2\Record Node 113\experiment1\recording1\events\Intan_RHD_USB-100.Rhythm Data\TTL' ;
csc_dir = '\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\test_data\SOM2_2025-10-22_15-21-16_opto_test_wide_probe2\Record Node 118'; 
phy_dir = '\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\Kilo_inter\SOM2_opto2_2';
vr_fname = []; %'/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D3/HF2b2_D3_2025-12-18_14-08-34.csv'; 
csc_idx = [3 11]; 
swr_ch = 1; 
save_name = 'SOM2_test_wide_probe2'; 
TTL = {'3', '7'};

%% HF2b2_H2_opto  [timing is off]
evts_dir = '\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\GoNoGo\HF2b2_H2\HF2b2_2026-01-12_09-38-52_H2_opto\Record Node 112\experiment1\recording1\events\Intan_RHD_USB-108.Rhythm Data\TTL' ;
csc_dir = '\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\GoNoGo\HF2b2_H2\HF2b2_2026-01-12_09-38-52_H2_opto\Record Node 117'; 
phy_dir = '\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\Kilo_inter\HF2b2_H2_opto';
vr_fname = []; %'/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D3/HF2b2_D3_2025-12-18_14-08-34.csv'; 
csc_idx = [3 11]; 
swr_ch = 1; 
save_name = 'HF2b2_H2_opto'; 
TTL = {'6' '7'};
%% convert for PC
if ispc
    usr =getenv('USERNAME'); 
    pc_name = ['C:\Users\' usr]; 
else
    usr =getenv('USER'); 
    pc_name = ['/Users/' usr]; 
end
    

evts_dir = regexprep([pc_name evts_dir], {'\', '/'}, {filesep, filesep});
csc_dir = regexprep([pc_name csc_dir], {'\', '/'}, {filesep, filesep});
phy_dir = regexprep([pc_name phy_dir], {'\', '/'}, {filesep, filesep});
if ~isempty(vr_fname)
vr_fname = regexprep([pc_name vr_fname], {'\', '/'}, {filesep, filesep});
end


%% plot the csc channels to check for good ones:
data = HF_preprocess(phy_dir, csc_dir, evts_dir, vr_fname, [1, 6, 11, 33:36, 48:52, 64:68]);

figure(1)
clf
hold on
offset =250;
for ii  = 1:size(data.csc.data, 1)
    plot(data.csc.tvec, data.csc.data(ii,:)+(offset*ii));

end
%% load the evts

data = HF_preprocess(phy_dir, csc_dir, evts_dir, vr_fname, csc_idx);

%% get the basic metrics

data = HF_metrics(data, 0, TTL); 

% remove cells with low firing rates
s_fr = [];
for ii = length(data.S_metrics):-1:1
    s_fr(ii) = data.S_metrics{ii}.fr;
end

data.S.t(s_fr<0.5) = [];
data.S.label(s_fr<0.5) = [];
data.S.usr(s_fr<0.5) = [];
data.S_metrics(s_fr<0.5) = [];


%% detect SWR

[data.swr.iv, data.swr.cfg] = MS_SWR_detector(data.csc, data.csc.label{swr_ch},1);
data.swr.chan = data.csc.label{swr_ch};
%% save the intermediate data
save(['/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Inter_data/' save_name '.mat'], "data")
%% figure showing all the events
c_ord = MS_linspecer(8);

% convert swr events to a rate.


figure(1010)
clf
hold on
plot(data.S);
plot(data.csc.tvec, (data.csc.data(swr_ch,:)/100)-2, 'b')
plot(data.licks.tvec, (data.licks.data./5)-10, 'Color', [0.6 .6 .6])

ylim([-10 length(data.S.t)+1])
for ii = 1:length(data.evts.label)
    if ~isempty(data.evts.t{ii})
        xline(data.evts.t{ii}(:),'Color', c_ord(ii,:), 'Label', data.evts.label{ii}, 'LabelColor',c_ord(ii,:))
    end
end


if isfield(data, 'vr')
    plot(data.vr.pos.tvec, data.vr.pos.data/5,'color', c_ord(3,:), 'LineWidth',3)
    xline(data.vr.evt.t{contains(data.vr.evt.label, 'Collision with Rwd1')}, 'c', 'LineStyle', '--')
end

if isfield(data, 'swr')
    swr_rate = MS_spike2rate(ts({IVcenters(data.swr.iv)}, {'swr'}), data.csc.tvec);
    yyaxis right
    plot(swr_rate.tvec, (swr_rate.data))
    ylim([0 25])
end


figure(1012)
clf
hold on
plot(data.S);
plot(data.csc.tvec, (data.csc.data(swr_ch,:)/100)-2, 'b')
plot(data.licks.tvec, (data.licks.data./5)-10, 'Color', [0.6 .6 .6])
ylim([-10 length(data.S.t)+1])
for ii = 1:length(data.evts.label)
    if ~isempty(data.evts.t{ii})
        xline(data.evts.t{ii}(:),'Color', c_ord(ii,:), 'Label', data.evts.label{ii}, 'LabelColor',c_ord(ii,:))
    end
end

if isfield(data, 'vr')
    plot(data.vr.pos.tvec, data.vr.pos.data/5, 'g', 'LineWidth',2)
    xline(data.vr.evt.t{contains(data.vr.evt.label, 'Collision with Rwd1')}, 'c', 'LineStyle', '--')
    xlim([data.vr.pos.tvec(1) data.vr.evt.t{contains(data.vr.evt.label, 'Collision with Rwd1')}(1)+5])
end

if isfield(data, 'swr')
    swr_rate = MS_spike2rate(ts({IVcenters(data.swr.iv)}, {'swr'}), data.csc.tvec);
    yyaxis right
    plot(swr_rate.tvec, (swr_rate.data))
    ylim([0 25])
end

%% peth per cell
c_red = [0.9153    0.2816    0.2878];


red_idx = find(contains(data.evts.label, TTL{1})); 
blue_idx = find(contains(data.evts.label, TTL{2})); 

window = [-.250 .250];
bin_s = .025;
evt_t = data.evts.t{red_idx} ;

% isolate events of a certain length
e_d = evt_t(2,:) - evt_t(1,:);

ITIs = unique(round(e_d, 3));

for iTi = 1:length(ITIs)

    evt_t = data.evts.t{red_idx} ;

    k_idx = round(e_d,3) == ITIs(iTi) ;

    evt_t(:,~k_idx) = [];
    %
    % cfg = [];
    % cfg.dt = 0.025;
    % cfg.gauss_window = .01;
    % cfg.gauss_sd = 0.0025;

    for iS = 1:length(data.S.t)

        % figure(iS+300)
        % clf

        S = data.S;
        S.t = [];
        S.t{1} = data.S.t{iS};
        S.label = [];
        S.label{1} =  data.S.label{iS};

        evt_t = sort(evt_t(1,:));

        cfg_peth = [];
        cfg_peth.window = [-.25 .25];
        cfg_peth. plot_type = 'raw';
        cfg_peth.dt = 0.001;
        cfg_peth.gauss_window = .025;
        cfg_peth.gauss_sd = .0025;
        cfg_peth.shuff = 500;
        cfg_peth.t_on = mode(e_d(k_idx));
        cfg_peth.rec_color = c_red;
        cfg_peth.plot = 'off';
        [peth_S{iS,iTi}, peth_IT{iS,iTi},peth_gau{iS,iTi}, ~, ~, ~, ~, ~, ~,peth_T{iS, iTi}] = SpikePETH_Shuff(cfg_peth, S, evt_t);
    end


end

%% plot the PETHS together
reds = hot(10);

for iS = 1:size(peth_gau,1)

    figure(iS);
    clf
    subplot(2,1,1)
    cla
    hold on;
    % Plot the PETH for each ITI
    offset = 0;
    for iT = 1:size(peth_gau,2)
        u_val = unique(peth_T{iS, iT});
        for iV = 1:length(u_val)
            this_idx = peth_T{iS, iT} == u_val(iV);
            if isempty(this_idx)
                continue
            end
            plot(peth_S{iS, iT}(this_idx), peth_T{iS, iT}(this_idx)+0.5 + offset,'.', 'color', reds(iT,:), 'MarkerSize', 10)
            % disp(mode(peth_T{iS, iT}(this_idx)+ offset))
        end
        y_lim = ylim;
        offset = y_lim(2);

    end
    ylabel('pulse #');


    subplot(2,1,2)
    hold on;
    % Plot the mean activity for each ITI
    for iT = 1:size(peth_gau,2)
        plot(peth_IT{iS,iT}, mean(peth_gau{iS,iT},2, 'omitnan')', 'Color', reds(iT,:), 'LineWidth', 1.5);

    end
    xlabel('Time (s)');
    ylabel('Firing rate (Hz)');
    % if data.

    title(['PETH for Cell ' num2str(iS) ' | Shank: ' num2str(data.S.usr{iS}.shank)...
        ' | x: ' num2str(round(data.S.usr{iS}.pos(1))), ' | y: ' num2str(round(data.S.usr{iS}.pos(2)))]);

    labels = [];
    for ii = 1:length(ITIs)
        labels{ii} = [num2str(ITIs(ii)*1000) ' ms | p = ' num2str(round(data.S_metrics{iS}.opto_red{ii}.pval, 3))];
    end

    legend(labels, 'Box', 'off')

end
red_labels = labels;

%% BLUE

clear peth_*
c_blue = [0.2878 0.2816 0.9153 ];

blue_idx = find(contains(data.evts.label, '7'));

window = [-.250 .250];
bin_s = .025;
evt_t = data.evts.t{blue_idx} ;

% isolate events of a certain length
e_d = evt_t(2,:) - evt_t(1,:);

ITIs = unique(round(e_d, 3));

for iTi = 1:length(ITIs)

    evt_t = data.evts.t{blue_idx} ;

    k_idx = round(e_d,3) == ITIs(iTi) ;

    evt_t(:,~k_idx) = [];

    for iS = 1:length(data.S.t)

        S = data.S;
        S.t = [];
        S.t{1} = data.S.t{iS};
        S.label = [];
        S.label{1} =  data.S.label{iS};

        evt_t = sort(evt_t(1,:));

        cfg_peth = [];
        cfg_peth.window = [-.25 .25];
        cfg_peth. plot_type = 'raw';
        cfg_peth.dt = 0.001;
        cfg_peth.gauss_window = .025;
        cfg_peth.gauss_sd = .0025;
        cfg_peth.shuff = 500;
        cfg_peth.t_on = mode(e_d(k_idx));
        cfg_peth.rec_color = c_blue;
        cfg_peth.plot = 'off';
        [peth_S{iS,iTi}, peth_IT{iS,iTi},peth_gau{iS,iTi}, ~, ~, ~, ~, ~, ~,peth_T{iS, iTi}] = SpikePETH_Shuff(cfg_peth, S, evt_t);
    end


end

%% BLUE plot the PETHS together
blues = cool(8);

for iS = 1:size(peth_gau,1)

    figure(iS+100);
    clf
    subplot(2,1,1)
    cla
    hold on;
    % Plot the PETH for each ITI
    offset = 0;
    for iT = 1:size(peth_gau,2)
        u_val = unique(peth_T{iS, iT});
        for iV = 1:length(u_val)
            this_idx = peth_T{iS, iT} == u_val(iV);
            if isempty(this_idx)
                continue
            end
            plot(peth_S{iS, iT}(this_idx), peth_T{iS, iT}(this_idx)+0.5 + offset,'.', 'color', blues(iT,:), 'MarkerSize', 10)
            % disp(mode(peth_T{iS, iT}(this_idx)+ offset))
        end
        y_lim = ylim;
        offset = y_lim(2);

    end
    ylabel('pulse #');


    subplot(2,1,2)
    hold on;
    % Plot the mean activity for each ITI
    for iT = 1:size(peth_gau,2)
        plot(peth_IT{iS,iT}, mean(peth_gau{iS,iT},2, 'omitnan')', 'Color', blues(iT,:), 'LineWidth', 1.5);

    end
    xlabel('Time (s)');
    ylabel('Firing rate (Hz)');
    % if data.

    title(['PETH for Cell ' num2str(iS) ' | Shank: ' num2str(data.S.usr{iS}.shank)...
        ' | x: ' num2str(round(data.S.usr{iS}.pos(1))), ' | y: ' num2str(round(data.S.usr{iS}.pos(2)))]);

    labels = [];
    for ii = 1:length(ITIs)
        labels{ii} = [num2str(ITIs(ii)*1000) ' ms | p = ' num2str(round(data.S_metrics{iS}.opto_blue{ii}.pval, 3))];
    end

    legend(labels, 'Box', 'off')

end

blue_labels = labels;
%% event triggered LFP average
figure(900)
clf

% RED
hold on
evt_t = data.evts.t{red_idx} ;

% isolate events of a certain length
e_d = evt_t(2,:) - evt_t(1,:);

ITIs = unique(round(e_d, 3));
win = data.csc.cfg.hdr{1}.SamplingFrequency / 4;

for iTi = 1:length(ITIs)

    evt_t = data.evts.t{red_idx} ;

    k_idx = round(e_d,3) == ITIs(iTi) ;

    evt_t(:,~k_idx) = [];

    evt_idx = nearest_idx3(evt_t, data.csc.tvec);

    sta_lfp{iTi} = [];
    for ii = length(evt_t):-1:1
        if (evt_idx(ii) + win) < length(data.csc.tvec)
            sta_lfp{iTi}(ii,:) = data.csc.data(swr_ch,evt_idx(ii) - win: evt_idx(ii)+win);
        end
    end
    plot(-.25:1/(data.csc.cfg.hdr{1}.SamplingFrequency):.25, mean(sta_lfp{iTi}), 'color', reds(iTi,:), 'LineWidth',1.5)

end

ylim([-40 40])
ylabel('Amplitude (\muV)')
set(gca, 'xtick', -.25:.05:.25)
legend(red_labels, 'Box', 'off')
xlabel('time from light (s)')
ylim([-100 100])


% BLUE
figure(901)
clf

hold on
evt_t = data.evts.t{blue_idx} ;

% isolate events of a certain length
e_d = evt_t(2,:) - evt_t(1,:);

ITIs = unique(round(e_d, 3));
win = data.csc.cfg.hdr{1}.SamplingFrequency / 4;

for iTi = 1:length(ITIs)

    evt_t = data.evts.t{red_idx} ;

    k_idx = round(e_d,3) == ITIs(iTi) ;

    evt_t(:,~k_idx) = [];

    evt_idx = nearest_idx3(evt_t, data.csc.tvec);

    sta_lfp{iTi} = [];
    for ii = length(evt_t):-1:1
        if (evt_idx(ii) + win) < length(data.csc.tvec)
            sta_lfp{iTi}(ii,:) = data.csc.data(swr_ch,evt_idx(ii) - win: evt_idx(ii)+win);
        end
    end
    plot(plot(-.25:1/(data.csc.cfg.hdr{1}.SamplingFrequency):.25, sta_lfp{iTi}, 'color', blues(iTi,:), 'LineWidth',.5))
    plot(-.25:1/(data.csc.cfg.hdr{1}.SamplingFrequency):.25, mean(sta_lfp{iTi}), 'color', blues(iTi,:), 'LineWidth',1.5)

end
ylim([-100 100])

ylabel('Amplitude (\muV)')
xlabel('time from light (s)')
set(gca, 'xtick', -.25:.05:.25)
legend(blue_labels, 'Box', 'off')


%% loop over sessions and collect the stats on the responsive cells

% grab the good intermediate sessions
int_fname= dir('/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Inter_data/*.mat');

keep_sess = {'HF2b2_D5_opto2.mat', 'HF2b2_D5.mat', 'HF2b2_D4'};
all_sess = [];
for ii = 1:length(int_fname)
    if contains(int_fname(ii).name, keep_sess)
        all_sess.(int_fname(ii).name(1:end-4)) = load([int_fname(ii).folder filesep int_fname(ii).name], 'data');
    end
end


% loop over cells and extract the responses
s_names = fieldnames(all_sess);
all_data = [];
for ii = length(s_names):-1:1

    this_data = all_sess.(s_names{ii}).data;
    cell = [];

    % prefill the responses
            cell.red_resp = NaN(length(this_data.S.t), 4); 
        cell.red_resp_dir = cell.red_resp; 
        cell.blue_resp = cell.red_resp;
        cell.blue_resp_dir = cell.red_resp; 

    for iC = length(this_data.S.t):-1:1

        cell.id(iC) = (ii*1000) + iC;
        cell.sess(iC) = ii;
        cell.fr(iC) = this_data.S_metrics{iC}.fr;
        cell.b_idx(iC) = this_data.S_metrics{iC}.burst_idx;
        cell.ISI(iC) = this_data.S_metrics{iC}.ISI;
        % cell.shank(iC) = this_data.S.usr{iC}.shank;
        cell.pos(iC,:) = this_data.S.usr{iC}.pos;


        for iT = length(this_data.S_metrics{iC}.opto_red):-1:1

            % if length(this_data.S_metrics{iC}.opto_red{iT}.pre) > 200% correct for cases where there are pulses with the same duration but different ITIs

            % red
            [~, cell.red_resp(iC, iT), ~, stats] = ttest(this_data.S_metrics{iC}.opto_red{iT}.pre, this_data.S_metrics{iC}.opto_red{iT}.post);
            cell.red_resp_dir(iC, iT) = sign(stats.tstat);
        end % end reds

        for iT = length(this_data.S_metrics{iC}.opto_blue):-1:1
            % blue
            [~, cell.blue_resp(iC, iT), ~, stats] = ttest(this_data.S_metrics{iC}.opto_blue{iT}.pre, this_data.S_metrics{iC}.opto_blue{iT}.post);
            cell.blue_resp_dir(iC, iT) = sign(stats.tstat);
        end % end blues
    end % loop over cells

    if isempty(all_data)
        all_data = cell;
    else
        % collect across sessions
        all_data.id = [all_data.id, cell.id];
        all_data.sess = [all_data.sess, cell.sess];
        all_data.fr = [all_data.fr,  cell.fr];
        all_data.b_idx = [all_data.b_idx, cell.b_idx];
        all_data.ISI = [all_data.ISI, cell.ISI];
        % all_data.shank = [all_data.shank, cell.shank];
        all_data.pos = [all_data.pos; cell.pos];

        all_data.red_resp = [all_data.red_resp; cell.red_resp];
        all_data.red_resp_dir = [all_data.red_resp_dir; cell.red_resp_dir];
        all_data.blue_resp = [all_data.blue_resp; cell.blue_resp];
        all_data.blue_resp_dir = [all_data.blue_resp_dir; cell.blue_resp_dir];
    end
end% session


%% make some plots
c_ord = MS_linspecer(5); 

figure(101)
clf

% K means 
subplot(1,2,1)
  [~, ~] = MS_kmean_scatter([all_data.fr', all_data.ISI', all_data.b_idx'], 3, [1,2,3], 150);
    xlabel('firing rate')
    ylabel('ISI')
    zlabel('bursting index')
    % axis square

    view([45 45 45])
% legend({'group 1', 'group 2', 'group 3'},'Location','northeast')
set(gca, 'XDir', 'reverse', 'yDir', 'reverse')

subplot(1,2,2)
cla
hold on
k_idx = logical(all_data.red_resp(:,2) < 0.05); 
% shanks = unique(all_data.shank);


% for ii = 1:length(shanks)
  
    this_idx = k_idx'; %& (all_data.shank == shanks(ii)); 
    scatter3(all_data.fr(this_idx), (all_data.ISI(this_idx)), all_data.b_idx(this_idx), 150,c_ord(ii,:), 'filled', 'Marker','o', 'DisplayName',' responsive')


    this_idx = ~k_idx'; %& (all_data.shank == shanks(ii)); 
    scatter3(all_data.fr(this_idx), (all_data.ISI(this_idx)), all_data.b_idx(this_idx), 150,c_ord(ii,:), 'filled', 'Marker','d', 'MarkerFaceAlpha',.25, 'DisplayName',['not responsive'])


% end
grid on
xlabel('Firing rate (Hz)')
ylabel('log ISI')
zlabel('Burst index')
legend('Location','northeast')
view([45 45 45])
set(gca, 'XDir', 'reverse', 'yDir', 'reverse')
% axis square

SetFigure([], gcf)
%%
figure(102)
cla
hold on

    scatter(all_data.pos(k_idx,1), all_data.pos(k_idx,2), 200,all_data.fr(k_idx), 'filled', 'Marker','o', 'DisplayName',['Shank ' num2str(ii) ' responsive'])
    scatter(all_data.pos(~k_idx,1), all_data.pos(~k_idx,2), 200,all_data.fr(~k_idx), 'filled','Marker','d', 'MarkerFaceAlpha',.25, 'DisplayName',['not responsive'])
cb = colorbar;
ylabel(cb,'Firing rate (Hz)','FontSize',16,'Rotation',90);
    ylim([-100 0])

    ylabel('Probe depth (microns)')
    xlabel('Probe width (microns)')
    text(250, 0, 'CA1', 'HorizontalAlignment','center', 'FontSize',22, 'VerticalAlignment','top')
text(0, 0, '<-   distal', 'HorizontalAlignment','left', 'FontSize',22, 'VerticalAlignment','top')
text(500, 0, 'intermediate   -> ', 'HorizontalAlignment','right', 'FontSize',22, 'VerticalAlignment','top')

SetFigure([], gcf)

%%
figure(103)
subplot(1,2,1)
r_10= (all_data.red_resp(:,1) < 0.05); 
r_50= (all_data.red_resp(:,2) < 0.05); 
r_100= (all_data.red_resp(:,3) < 0.05); 
r_250= (all_data.red_resp(:,4) < 0.05); 
n_cells = length(all_data.red_resp(:,4)); 

ax = piechart([sum(r_10 | r_50 | r_100 | r_250)/n_cells, ...        
    sum(r_10 & ~r_50 &~r_100 &~r_250)/n_cells,...
    sum(~r_10 & r_50 &~r_100 &~r_250)/n_cells,...
    sum(~r_10 & ~r_50 & r_100 &~r_250)/n_cells,...
    sum(~r_10 & ~r_50 &~r_100 & r_250)/n_cells,...
    sum(~r_10 & ~r_50 &~r_100 &~r_250)/n_cells,...
    ]*100, {'any', '10ms', '50ms', '100ms', '250ms', 'non-repsonive'}, 'fontsize', 12);

red_map = hot(12); 
colororder([red_map(1:2:10,:); .7 .7 .7])

title('Red light responses')

% 
% subplot(1,2,2)
% r_10= (all_data.blue_resp(:,1) < 0.05); 
% r_50= (all_data.blue_resp(:,2) < 0.05); 
% r_100= (all_data.blue_resp(:,3) < 0.05); 
% r_250= (all_data.blue_resp(:,4) < 0.05); 
% n_cells = length(all_data.blue_resp(:,4)); 
% 
% ax = piechart([sum(r_10 | r_50 | r_100 | r_250)/n_cells, ...
%     sum(r_10 & ~r_50 &~r_100 &~r_250)/n_cells,...
%     sum(~r_10 & r_50 &~r_100 &~r_250)/n_cells,...
%     sum(~r_10 & ~r_50 & r_100 &~r_250)/n_cells,...
%     sum(~r_10 & ~r_50 &~r_100 & r_250)/n_cells,...
%     sum(~r_10 & ~r_50 &~r_100 &~r_250)/n_cells,...
%     ]*100, {'any', '10ms', '50ms', '100ms', '250ms', 'non-repsonive'}, 'fontsize', 12);
% 
% title('Blue light responses')
% ax.Colormap = ([MS_linspecer(5); .7 .7 .7]);
% colororder([MS_linspecer(5); .7 .7 .7])