%% sandbox_HF_OE_Peth


%% HF2b2_D3
evts_dir = ('/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D5/HF2b2_2025-12-18_13-47-24_D3_2_opto_only/Record Node 112/experiment1/recording1/events/Intan_RHD_USB-108.Rhythm Data/TTL') ;
csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D3/HF2b2_2025-12-18_13-47-24_D3_2_opto_only/Record Node 117'; 
phy_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Kilo_inter/HF2b2_D3';
vr_fname = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D3/HF2b2_D3_2025-12-18_14-08-34.csv'; 


%% HF2b2_D5
evts_dir = ('/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D5/HF2b2_2026-01-02_13-50-31_D5/Record Node 112/experiment1/recording1/events/Intan_RHD_USB-108.Rhythm Data/TTL') ;
csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D5/HF2b2_2026-01-02_13-50-31_D5/Record Node 117'; 
phy_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Kilo_inter/HF2b2_D5';
vr_fname = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/GoNoGo/HF2b2_D5/h2b2_VR_D5_2026-01-02_14-09-33.csv'; 
csc_idx = [1 6 11]; 
swr_ch = 2; 

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

data = HF_metrics(data); 

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

red_idx = find(contains(data.evts.label, '6')); 
blue_idx = find(contains(data.evts.label, '7')); 

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
    plot(-.25:1/(data.csc.cfg.hdr{1}.SamplingFrequency):.25, mean(sta_lfp{iTi}), 'color', blues(iTi,:), 'LineWidth',1.5)

end
ylim([-40 40])

ylabel('Amplitude (\muV)')
xlabel('time from light (s)')
set(gca, 'xtick', -.25:.05:.25)
legend(blue_labels, 'Box', 'off')


%% loop over sessions and collect the stats on the responsive cells

% grab the good intermediate sessions
int_fname= dir('/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Inter_data/*.mat'); 

keep_sess = {'HF2b2_D5_opto2.mat', 'HF2b2_D5.mat'}; 

for ii = 1:length(int_fname)
    if contains(keep_sess, int_fname(ii).name)
        all_data.(int_fname(ii).name(1:end-4)) = load([int_fname(ii).folder filesep int_fname(ii).name]);
    end
end


% loop over cells and extract the responses 
s_names = fieldnames(all_data); 

for ii = length(s_names):-1:1

    this_sess = all_data.(s_names{ii}); 

    for iC = length(this_sess.S.t):-1:1
        
        for iO = 1:2  % loop over 



    end % loop over cells
end% session
