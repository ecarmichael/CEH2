function [outputS, outputIT, outputGau,outputGau_shuf, pre_stim_means, post_stim_means, pre_stim_std, post_stim_std, z_vals, outputT] = SpikePETH_Shuff(cfg_in, S,t,varargin)
%% SpikePETH_Shuff: computes the perievent histogram for spike data "S" at events
%             "t".  Outputs
%
%
%
%          Inputs:
%           - cfg_in [struct]: contains configuration paramters
%           - S      [TS]      Spike timestamp data
%           - t      [n x T]   timestamps for events
%          Outputs:
%           -
%           -
%           -
% based on spikePETH by MvdM
% modified by EC to match mvdmlab codebase- 2017-05-01

%% set defaults
cfg_def.shuff = 5000; % give a value for shuffle. will skip if empty
cfg_def.window = [-2 5];
cfg_def.dt = 0.001;
cfg_def.t_on = [];
cfg_def.excessBounds = 1;
cfg_def.outputGrid = 0;
cfg_def.evt_color_mat = repmat([0 0 0], length(t),1);
cfg_def.rec_color = [4,172,218]./255;
cfg_def.linewidth = 2;
cfg_def.markersize = 5;
cfg_def.color = [0.3639    0.5755    0.7484];
cfg_def.binsize = cfg_def.dt; % used for gaussian kernal.  select a small bin size for good time resolution
cfg_def.waves = [];
cfg_def.contrast_waves = [];
cfg_def.c_ord = linspecer(4);
cfg_def.gauss_window = 1;
cfg_def.gauss_sd = 0.02;
cfg_def.plot = 'on'; % turn output 'on' or 'off';
cfg_def.plot_type = 'zscore'; % could also be zscore.
cfg_def.z_mean = [];
cfg_def.z_std = [];
cfg = ProcessConfig2(cfg_def, cfg_in);

cfg.binsize = cfg.dt;
extract_varargin;

c_ord = linspecer(5);

%% compute the spike PETH
nT = length(t);

outputS = [];
outputT = [];
outputGau = [];
% outputID = repmat(inf, nT, diff(cfg.window)/cfg.dt+1);
outputIT = linspace(cfg.window(1), cfg.window(2), diff(cfg.window)/cfg.dt+1);

if cfg.outputGrid
    xbin = cfg.window(1):cfg.dt:cfg.window(2);
    outputG = zeros(nT,length(xbin)-1);
end

% set up gau kernal
gauss_window = cfg.gauss_window./cfg.binsize; % 1 second window
gauss_SD = cfg.gauss_sd./cfg.binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./cfg.binsize; % normalize by binsize

% convolve with gaussian for firing rate

% plot(tbin_centers,S_gau_sdf,'g');

for iT = nT:-1:1
    S0 = restrict(S, t(iT)+cfg.window(1)-cfg.excessBounds, t(iT)+cfg.window(2)+cfg.excessBounds);
    if length(S0.t{1}) > 0
        S0 = restrict(S0, t(iT)+cfg.window(1), t(iT)+cfg.window(2));

        outputT = [outputT; repmat(iT, length(S0.t{1}),1)];
        outputS = [outputS; S0.t{1}-t(iT)];

        %convolve with gaussian for firing rate.
        tbin_edges = t(iT)+cfg.window(1):cfg.binsize:t(iT)+cfg.window(2);
        tbin_centers = tbin_edges(1:end-1)+cfg.binsize/2;
        spk_count = histc(S0.t{1},tbin_edges);
        spk_count = spk_count(1:end-1);


        S_gau_sdf = conv2(spk_count,gk,'same'); % convolve with gaussian window
        if size(S_gau_sdf,1) >1
            S_gau_sdf = S_gau_sdf';
        end
        outputGau(:,iT) = S_gau_sdf;


        if cfg.outputGrid

            temp = histc(S0.t{1}-t(iT),xbin); temp = temp(1:end-1);
            if ~isempty(temp)
                outputG(iT,:) = temp;
            end

        end

    end
end
outputIT = outputIT(1:end-1);


%% check if there are any spike\
if isempty(outputT)
    z_vals = nan(size(outputIT))';
    outputS = nan(size(outputIT))';
    outputGau = nan(size(outputIT))';
    outputS_shuf = [];
    outputT_shuf = [];
    outputGau_shuf = [];
    mean_S_gau = nan(size(outputIT))';
    pre_stim_means = nan(size(t));
    post_stim_means = nan(size(t));
    pre_stim_std = nan(size(t));
    post_stim_std =nan(size(t));
    disp('No spikes: filling with NaNs')
    return
end

%% shuffles

rng(101, 'twister'); % for reproducibility.
shuff_t =  ((max(S.t{1})-cfg.window(2))-(min(S.t{1}) + abs(cfg.window(1)))).*rand(cfg.shuff,1) + (min(S.t{1}) + abs(cfg.window(1))); %random time points between the first and last spike with window size.
outputS_shuf = [];
outputT_shuf = [];
outputGau_shuf = [];
% outputID = repmat(inf, nT, diff(cfg.window)/cfg.dt+1);

if cfg.outputGrid
    xbin = cfg.window(1):cfg.dt:cfg.window(2);
    outputG = zeros(nT,length(xbin)-1);
end



for iShuf = cfg.shuff:-1:1

    S0 = restrict(S, shuff_t(iShuf)+cfg.window(1)-cfg.excessBounds, shuff_t(iShuf)+cfg.window(2)+cfg.excessBounds);
    if length(S0.t{1}) > 0
        S0 = restrict(S0, shuff_t(iShuf)+cfg.window(1), shuff_t(iShuf)+cfg.window(2));


        outputT_shuf = [ repmat(iShuf, length(S0.t{1}),1);  outputT_shuf];
        outputS_shuf = [ S0.t{1}-shuff_t(iShuf); outputS_shuf];

        %convolve with gaussian for firing rate.
        tbin_edges = shuff_t(iShuf)+cfg.window(1):cfg.binsize:shuff_t(iShuf)+cfg.window(2);
        tbin_centers = tbin_edges(1:end-1)+cfg.binsize/2;
        spk_count = histc(S0.t{1},tbin_edges);
        spk_count = spk_count(1:end-1);

        gauss_window = cfg.gauss_window./cfg.binsize; % 1 second window
        gauss_SD = cfg.gauss_sd./cfg.binsize; % 0.02 seconds (20ms) SD
        gk = gausskernel(gauss_window,gauss_SD); gk = gk./cfg.binsize; % normalize by binsize
        S_gau_sdf = conv2(spk_count,gk,'same'); % convolve with gaussian window
        if size(S_gau_sdf,1) >1
            S_gau_sdf = S_gau_sdf';
        end
        outputGau_shuf(:,iShuf) =  S_gau_sdf;


        if cfg.outputGrid

            temp = histc(S0.t{1}-t(iT),xbin); temp = temp(1:end-1);
            if ~isempty(temp)
                outputG(iT,:) = temp;
            end

        end

    end % end check for any spikes
end % end shuffles.


%% get the zscore and get means;

mean_S_gau = nanmean(outputGau,2); % get the mean gaussian smoothed firing rate

mean_S_gau_z = nanmean(outputGau,2); % get the mean gaussian smoothed firing rate
mean_S_gau_z = (mean_S_gau_z - nanmean(outputGau_shuf,2))./nanstd(outputGau_shuf, [], 2);

idx = nearest_idx(0, outputIT); % get the event time index

pre_stim_means = nanmean(outputGau(1:idx-1,:),1);  % get the mean of the gau smoothed firing rate before the event.
post_stim_means = nanmean(outputGau(idx:end,:),1); % % get the mean of the gau smoothed FR after the event.

pre_stim_std = nanstd(outputGau(1:idx-1,:),[],1);  % get the mean of the gau smoothed firing rate before the event.
post_stim_std = nanstd(outputGau(idx:end,:),[],1); % % get the mean of the gau smoothed FR after the event.

z_vals = (mean_S_gau - mean(mean_S_gau(1:idx)))./mean(mean_S_gau(1:idx));




%% display
clf
if  strcmp(cfg.plot, 'on')


    % spike raster
    subplot(2,1,1);
    % 	imagesc(window,[1 nT], outputID);
    % 	colormap(1-0.25*gray);
    hold on;
    u_val = unique(outputT);
    for iV = 1:length(u_val)
        this_idx = outputT == u_val(iV);
        if isempty(this_idx)
            continue
        end
        plot(outputS(this_idx), outputT(this_idx)+0.5,'.', 'color', cfg.evt_color_mat(u_val(iV),:), 'MarkerSize', cfg.markersize)
    end

    %     plot(outputS, outputT+0.5, '.k', 'MarkerSize', 5);
    xlabel('peri-event (sec)');
    ylabel('Event #');
    if nT < 2
        ylim([0 1])
    else
        ylim([1 nT])
    end
    xlim(cfg.window);
    hold on
    if ~isempty(cfg.t_on)
        rectangle('position', [0 1 cfg.t_on  nT], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
    elseif (size(t,2) > 1) || (size(t,1) == 1)
        rectangle('position', [0 1 0.001  nT], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
    else
        rectangle('position', [0 1 abs(mode(t(:,2)-t(:,1)))  nT], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
    end

    %Reverse the stacking order so that the patch overlays the line
    set(gca, 'Children',flipud(get(gca, 'Children')))
    %% add in the wave forms
    if ~isempty(cfg.waves)
        axes('Position', [.8 .85 0.1 .075])
        for iV = 1:4
            hold on
            plot(cfg.waves.xrange(:,iV), cfg.waves.mWV(:,iV), 'color', cfg.c_ord(iV, :), 'linewidth', 2)
            set(gca, 'visible', 'off')
        end
    end

    %%
    % bar graph
    % subplot(3,2,3);
    % m = histc(outputS, outputIT);
    % bar(outputIT,m/cfg.dt/length(t));
    % % 	x = outputIT;
    % % 	m = nanmean(1./outputID);
    % % 	se =  nanstd(1./outputID)/sqrt(nT+1);
    % % 	plot(x,m,'b',x,m+se,'r:',x,m-se,'r:');
    % set(gca, 'XLim', cfg.window);
    % ylabel('FR (Hz)')
    % xlabel('peri-event (sec)');

    % mean frequency line
    subplot(2,1,2);
    yyaxis left
    hold on
    %     outputITG= linspace(cfg.window(1), cfg.window(2), diff(cfg.window)/gauss_window+1);

    % se_S_gau = nanstd(outputGau,2)/sqrt(nT+1);
    % plot(outputIT(1:end-1),mean_S_gau, 'b',outputIT(1:end-1),mean_S_gau+se_S_gau, 'b:',outputIT(1:end-1),mean_S_gau-se_S_gau, 'b:' )


    if strcmp(cfg.plot_type, 'zscore')
        plot(outputIT, z_vals,'color', 'k', 'linewidth', cfg.linewidth)
        ylabel('Pre-event zscore')
        ylim([min(z_vals) max(z_vals)]);
        if ~isempty(cfg.shuff)
            plot(outputIT,nanmean(zscore(outputGau_shuf,[], 'all'),2), '--', 'color', [0.3 .3 .3])
        end
        if ~(max(mean_S_gau)) ==0
            if ~isempty(cfg.t_on)
                rectangle('position', [0 min(z_vals) cfg.t_on  (max(z_vals) - min(z_vals))*10], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
            elseif (size(t,2) > 1) || (size(t,1) == 1)
                rectangle('position', [0 min(z_vals) 0.001  (max(z_vals) - min(z_vals))*10], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
            else
                rectangle('position', [0 min(z_vals) abs(mode(t(:,2)-t(:,1)))  (max(z_vals) - min(z_vals))*10], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
            end
        end
        set(gca, 'Children',flipud(get(gca, 'Children')))

        x_lims = xlim;
        y_lims = ylim;

        mean_gau =nanmean(z_vals,2);

        text(x_lims(1), y_lims(2)*.9, ['Pre mean: ' num2str(mean(mean_gau(1:idx-1)), 2) '+/-' num2str(std(mean_gau(1:idx-1)),2) 'Hz'], 'fontweight', 'bold', 'fontsize', 12, 'color',c_ord(1,:) )
        text(x_lims(1), y_lims(2)*.7, ['Post mean: ' num2str(mean(mean_gau(idx:end)), 2)  '+/-' num2str(std(mean_gau(idx:end)),2) 'Hz' ], 'fontweight', 'bold', 'fontsize', 12, 'color',c_ord(2,:))



    elseif strcmp(cfg.plot_type, 'shuff_zscore')
        plot(outputIT, mean_S_gau_z,'color', 'k', 'linewidth', cfg.linewidth)
        ylabel('Pre-event shuff zscore')
        ylim([min(mean_S_gau_z) max(mean_S_gau_z)]);
        %         if ~isempty(cfg.shuff)
        %             plot(outputIT,nanmean(zscore(outputGau_shuf,[], 'all'),2), '--', 'color', [0.3 .3 .3])
        %         end
        if ~(max(mean_S_gau_z)) ==0
            if ~isempty(cfg.t_on)
                rectangle('position', [0 min(mean_S_gau_z) cfg.t_on  (max(mean_S_gau_z) - min(mean_S_gau_z))*10], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
            elseif (size(t,2) > 1) || (size(t,1) == 1)
                rectangle('position', [0 min(mean_S_gau_z) 0.001  (max(mean_S_gau_z) - min(mean_S_gau_z))*10], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
            else
                rectangle('position', [0 min(mean_S_gau_z) abs(mode(t(:,2)-t(:,1)))  (max(mean_S_gau_z) - min(mean_S_gau_z))*10], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
            end
        end
        x_lims = xlim;
        y_lims = ylim;
        set(gca, 'Children',flipud(get(gca, 'Children')))

        mean_gau =nanmean(mean_S_gau_z,2);

        text(x_lims(1), y_lims(2)*.9, ['Pre mean: ' num2str(mean(mean_gau(1:idx-1)), 2) '+/-' num2str(std(mean_gau(1:idx-1)),2) 'Hz'], 'fontweight', 'bold', 'fontsize', 12, 'color',c_ord(1,:) )
        text(x_lims(1), y_lims(2)*.7, ['Post mean: ' num2str(mean(mean_gau(idx:end)), 2)  '+/-' num2str(std(mean_gau(idx:end)),2) 'Hz' ], 'fontweight', 'bold', 'fontsize', 12, 'color',c_ord(2,:))


    else
        plot(outputIT, mean_S_gau,'color', 'k', 'linewidth', cfg.linewidth)
        xlim(cfg.window);
        ylabel('firing rate (Hz)');
        if ~isempty(cfg.shuff)
            plot(outputIT,nanmean(outputGau_shuf,2), '--', 'color', [0.3 .3 .3])
        end

        if ~(max(mean_S_gau)) ==0
            ylim([min(mean_S_gau) max(mean_S_gau)])
            if ~isempty(cfg.t_on)
                rectangle('position', [0 min(mean_S_gau) cfg.t_on  (max(mean_S_gau) - min(mean_S_gau))*10], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
            elseif (size(t,2) > 1) || (size(t,1) == 1)
                rectangle('position', [0 min(mean_S_gau) 0.001  abs(max(mean_S_gau))*10], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
            else
                rectangle('position', [0 min(mean_S_gau) abs(mode(t(:,2)-t(:,1)))  abs(max(mean_S_gau))*10], 'facecolor', [cfg.rec_color 0.5], 'edgecolor', [cfg.rec_color 0.5])
            end
        end
        set(gca, 'Children',flipud(get(gca, 'Children')))


        % add in pre and post event means
        x_lims = xlim;
        y_lims = ylim;

        mean_gau =nanmean(outputGau,2);

        text(x_lims(1), y_lims(2)*.9, ['Pre mean: ' num2str(mean(mean_gau(1:idx-1)), 2) '+/-' num2str(std(mean_gau(1:idx-1)),2) 'Hz'], 'fontweight', 'bold', 'fontsize', 12, 'color',c_ord(1,:) )
        text(x_lims(1), y_lims(2)*.7, ['Post mean: ' num2str(mean(mean_gau(idx:end)), 2)  '+/-' num2str(std(mean_gau(idx:end)),2) 'Hz' ], 'fontweight', 'bold', 'fontsize', 12, 'color',c_ord(2,:))

    end % end plot type

end % end plotting


if ~exist('z_vals')
    z_vals = nan(size(outputIT))';
end

if ~exist('outputGau')
    outputGau = nan(size(outputIT))';
end

if ~exist('outputGau_shuf')
    outputGau_shuf = nan(size(outputIT))';
end

if ~exist('pre_stim_means')
    pre_stim_means = nan(size(t));
end

if ~exist('post_stim_means')
    post_stim_means = nan(size(t));
end

if ~exist('pre_stim_std')
    pre_stim_std = nan(size(t));
end

if ~exist('post_stim_std')
    post_stim_std = nan(size(t));
end

end % end function
