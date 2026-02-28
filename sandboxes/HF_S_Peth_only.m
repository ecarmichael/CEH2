%% PETH only sandbox

%% JAWS_D1_1
evts_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/Jaws_test/JAWS22026-02-25_14-09-12_D1/Record Node 101/experiment1/recording1/events/Acquisition_Board-100.acquisition_board/TTL';
csc_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/Jaws_test/JAWS22026-02-25_14-09-12_D1/Record Node 117'; 
phy_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Kilo_inter/Jaws2_D1_1/kilosort4';
vr_fname = ''; 
csc_idx = [13]; 
swr_ch = 1; 
save_name = 'HF3b2_TFC_D4'; 
TTL = {'6'};

%% JAWS_D1_2
evts_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/Jaws_test/JAWS22026-02-25_14-21-05_D1_1/Record Node 101/experiment1/recording1/events/Acquisition_Board-100.acquisition_board/TTL';
csc_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/Jaws_test/JAWS22026-02-25_14-21-05_D1_1/Record Node 117'; 
phy_dir = '/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Kilo_inter/Jaws2_D1_2/kilosort4';
vr_fname = ''; 
csc_idx = [13]; 
swr_ch = 1; 
save_name = 'HF3b2_TFC_D4'; 
TTL = {'6'};

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


%%


params = OE_load_params(phy_dir);

data.S = OE_phy2TS(phy_dir, params);

ts_prime = readNPY([phy_dir filesep 'timestamps.npy']); 

data.evts = OE_load_binary_evts(evts_dir, ts_prime(1));


% make a colour map for the shanks
s_pos = []; 
for ii = length(data.S.usr):-1:1; s_pos(ii,:) = data.S.usr{ii}.pos; end

[snk_u,~,  s_idx] = unique(round(s_pos(:,1)/100)); % group by 100s of micros

c_ord = MS_linspecer(8);
spkColor = ones(length(data.S.t),3); 
for ii = 1:length(s_idx)
    spkColor(ii,:) = c_ord(s_idx(ii),:); 
end

[~, sort_idx] = sort(s_idx); 

data.S.t = data.S.t(sort_idx); 
data.S.label = data.S.label(sort_idx); 
data.S.usr = data.S.usr(sort_idx); 

spkColor = spkColor(sort_idx,:); 

%% quick metrics
S_metrics = []; 
for iS = length(data.S.t):-1:1
    S_metrics.fr(iS) = length(data.S.t{iS})./(ts_prime(end)-ts_prime(1));
    S_metrics.ISI(iS) = mean(diff(data.S.t{iS}));
    S_metrics.burst_idx(iS) = sum(diff(data.S.t{iS}) < .010) / sum(diff(data.S.t{iS}) > .010);
end


data_in = [S_metrics.fr; S_metrics.ISI; S_metrics.burst_idx]'; 
figure(1)
clf
[g_idx] = MS_kmean_scatter(data_in, 3, 1:length(data_in), 100); 

% scatter3(S_metrics.fr, S_metrics.ISI, S_metrics.burst_idx, 150, spkColor, 'filled')
xlabel('firing rate (Hz)')
ylabel('ISI (log)')
zlabel('burst index')
set(gca, 'YScale', 'log')

%% figure showing all the events

figure(1010)
clf
hold on

cfg = []; cfg.openNewFig = 0;
cfg.spkColor = spkColor; 
h = MultiRaster(cfg,data.S);

ylim([0 length(data.S.t)+1])
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
cfg = []; cfg.openNewFig = 0;
cfg.spkColor = spkColor; 

h = MultiRaster(cfg,data.S);
ylim([ 0 length(data.S.t)+1])
xlim([data.evts.t{1}(1) data.evts.t{1}(1)+30])
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

window = [-1 1];
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
        cfg_peth.window = window;
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
reds = hot(8);
reds = reds(2:end,:); 

for iS = 1:size(peth_gau,1)

    figure(iS);
    clf
    subplot(3,1,1:2)
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

            % % line plot
            % nSpikes = length(peth_T{iS, iT}(this_idx));
            % xvals = peth_T{iS, iT}(this_idx)';
            % yvals = iC.*(1,nSpikes);
            % 
            % xvals = xvals(:);
            % yvals = yvals(:);
            % 
            % if isempty(xvals) && isempty(yvals)
            %     %             h(iC) = nan;
            %     h{iC} = nan;
            % else
            %     %             h(iC) = plot(xvals,yvals,'Color',cfg.spkColor(iC,:),'LineWidth',cfg.LineWidth);
            %     plot(xvals,yvals,'.','Color',cfg.spkColor(iC,:));
            %     h{iC} = get(gca);
            % end
            % plot([peth_S{iS, iT}(this_idx)*1000, peth_S{iS, iT}(this_idx)*1000]', [peth_T{iS, iT}(this_idx)-0.5 + iT,peth_T{iS, iT}(this_idx)+0.5 + iT]', 'color', reds(iT,:), 'LineWidth',4)

            plot(peth_S{iS, iT}(this_idx)*1000, peth_T{iS, iT}(this_idx)+0.5 + offset,'.', 'Color', reds(iT,:), 'MarkerSize', 10)
            % disp(mode(peth_T{iS, iT}(this_idx)+ offset))
        end
        y_lim = ylim;
        offset = y_lim(2);

    end
    set(gca, 'XTicklabel', [])
    ylim([0 length(u_val)+1])
    ylabel('pulse #');
    xlim([-100 250])


    subplot(3,1,3)
    hold on;
       labels = {};
    for ii = 1:length(ITIs)
        rectangle('Position',[0, 0,ITIs(ii)*1000,  max(mean(peth_gau{iS,iT},2, 'omitnan')')*1.1], 'FaceColor',reds(iT,:), 'FaceAlpha',.2, 'EdgeColor','none')
        labels{ii} = [num2str(ITIs(ii)*1000) ' ms']; % | p = ' num2str(round(data.S_metrics{iS}.opto_red{ii}.pval, 3))];
    end

    % Plot the mean activity for each ITI
    for iT = 1:size(peth_gau,2)
        plot(peth_IT{iS,iT}*1000, mean(peth_gau{iS,iT},2, 'omitnan')', 'Color', reds(iT,:), 'LineWidth', 1.5);

    end
    xlabel('Time ms)');
    ylabel('Firing rate (Hz)');

    % if data.

    title(['PETH for Cell ' num2str(iS)])% ' | Shank: ' num2str(data.S.usr{iS}.shank)...
        %' | x: ' num2str(round(data.S.usr{iS}.pos(1))), ' | y: ' num2str(round(data.S.usr{iS}.pos(2)))]);

 
    legend(labels, 'Box', 'off')
    xlim([-100 250])
    ylim([0 max(mean(peth_gau{iS,iT},2, 'omitnan')')*1.1])


    % % same thing but zoomned in
    %     subplot(4,1,3)
    % cla
    % hold on;
    % % Plot the PETH for each ITI
    % offset = 0;
    % for iT = 1:size(peth_gau,2)
    %     u_val = unique(peth_T{iS, iT});
    %     for iV = 1:length(u_val)
    %         this_idx = peth_T{iS, iT} == u_val(iV);
    %         if isempty(this_idx)
    %             continue
    %         end
    %         plot(peth_S{iS, iT}(this_idx)*1000, peth_T{iS, iT}(this_idx)+0.5 + offset,'.', 'color', reds(iT,:), 'MarkerSize', 10)
    %         % disp(mode(peth_T{iS, iT}(this_idx)+ offset))
    %     end
    %     y_lim = ylim;
    %     offset = y_lim(2);
    % 
    % end
    % ylabel('pulse #');
    % xlim([-150 150])
    % 
    % subplot(4,1,4)
    % hold on;
    % % Plot the mean activity for each ITI
    % for iT = 1:size(peth_gau,2)
    %     plot(peth_IT{iS,iT}*1000, mean(peth_gau{iS,iT},2, 'omitnan')', 'Color', reds(iT,:), 'LineWidth', 1.5);
    % 
    % end
    % xlabel('Time ms)');
    % ylabel('Firing rate (Hz)');
    % % if data.
    % 
    % title(['PETH for Cell ' num2str(iS)])% ' | Shank: ' num2str(data.S.usr{iS}.shank)...
    %     %' | x: ' num2str(round(data.S.usr{iS}.pos(1))), ' | y: ' num2str(round(data.S.usr{iS}.pos(2)))]);
    % 
    % labels = {};
    % for ii = 1:length(ITIs)
    %     labels{ii} = [num2str(ITIs(ii)*1000) ' ms']; % | p = ' num2str(round(data.S_metrics{iS}.opto_red{ii}.pval, 3))];
    % end
    % xlim([-150 150])
    % 
    % legend(labels, 'Box', 'off')

end
red_labels = labels;


%% get teh xcorr

xbin_centers = -10-0.01:0.01:10+0.01; % first and last bins are to be deleted later
%ac = zeros(size(xbin_centers));
ac = cell(1,length(data.S.t)); 

parfor iC = 1:length(data.S.t)
    fprintf('%f.0,', iC)
    ac{iC} = zeros(size(xbin_centers));

    for iSpk = 1:length(data.S.t{iC})

        relative_spk_t = data.S.t{iC} - data.S.t{iC}(iSpk);

        ac{iC} = ac{iC} + hist(relative_spk_t,xbin_centers); % note that hist() puts all spikes outside the bin centers in the first and last bins! delete later.

    end
ac{iC} = ac{iC}(2:end-1);
ac{iC}(zero_idx) = 0;
end
fprintf('\n')

xbin = xbin_centers(2:end-1); % remove unwanted bins
zero_idx = find(xbin == 0);


%% plot the autocorr

ac_ord = MS_linspecer(length(ac)); 
figure(201)
clf
hold on
for ii = 1:length(ac)


    this_data = ac{ii};
    z_idx = nearest_idx(0, xbin); 
    this_data(z_idx - 2:z_idx+2) = NaN; 
    plot(xbin, (this_data./max(this_data))+ii, 'color', ac_ord(ii,:), 'LineWidth',2)

end
