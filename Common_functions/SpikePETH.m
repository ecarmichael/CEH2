function [outputS, outputT, outputGau, pre_stim_mean, post_stim_mean] = SpikePETH(cfg_in, S,t,varargin)
%% SpikePETH: computes the perievent histogram for spike data "S" at events
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
cfg_def.window = [-2 5];
cfg_def.dt = 0.00025;
cfg_def.excessBounds = 1;
cfg_def.outputGrid = 0;
cfg_def.evt_color = [4,172,218]./255;
cfg_def.binsize = cfg_def.dt; % used for gaussian kernal.  select a small bin size for good time resolution
cfg_def.waves = [];
cfg_def.contrast_waves = [];
cfg_def.c_ord = linspecer(4);
cfg_def.gauss_window = .001; 
cfg_def.gauss_sd = 0.0002;
cfg = ProcessConfig2(cfg_def, cfg_in);

extract_varargin;

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

% convolve with gaussian for firing rate

% plot(tbin_centers,S_gau_sdf,'g');

for iT = 1:nT
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
        
        gauss_window = cfg.gauss_window./cfg.binsize; % 1 second window
        gauss_SD = cfg.gauss_sd./cfg.binsize; % 0.02 seconds (20ms) SD
        gk = gausskernel(gauss_window,gauss_SD); gk = gk./cfg.binsize; % normalize by binsize
        S_gau_sdf = conv2(spk_count,gk,'same'); % convolve with gaussian window
        if size(S_gau_sdf,1) >1
            S_gau_sdf = S_gau_sdf';
        end
        outputGau = [outputGau; S_gau_sdf];
        
        
        if cfg.outputGrid
            
            temp = histc(S0.t{1}-t(iT),xbin); temp = temp(1:end-1);
            if ~isempty(temp)
                outputG(iT,:) = temp;
            end
            
        end
        
    end
end
%% check if there are any spikes
if isempty(outputT)
    disp('No spikes')
    return
end
%% display
clf

% spike raster
subplot(2,1,1);
% 	imagesc(window,[1 nT], outputID);
% 	colormap(1-0.25*gray);
% 	hold on;
plot(outputS, outputT+0.5, 'k.', 'MarkerSize', 5);
xlabel('peri-event (sec)');
ylabel('Event #');
ylim([1 nT])
xlim(cfg.window);
hold on
if size(t,2) > 1
rectangle('position', [0 1 abs(mode(t(:,2)-t(:,1)))  nT], 'facecolor', [cfg.evt_color 0.5], 'edgecolor', [cfg.evt_color 0.5])
else
rectangle('position', [0 1 0.001  nT], 'facecolor', [cfg.evt_color 0.5], 'edgecolor', [cfg.evt_color 0.5])
end
%% add in the wave forms
if ~isempty(cfg.waves)
    for ii = 1:4
        axes('Position', [(.65+(.05*ii)) .8 0.05 .1])
        plot(cfg.waves.mWV(:,ii), 'color', cfg.c_ord(ii, :))
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
mean_S_gau = nanmean(outputGau,1);
% se_S_gau = nanstd(outputGau,2)/sqrt(nT+1);
% plot(outputIT(1:end-1),mean_S_gau, 'b',outputIT(1:end-1),mean_S_gau+se_S_gau, 'b:',outputIT(1:end-1),mean_S_gau-se_S_gau, 'b:' )
plot(outputIT(1:end-1), mean_S_gau)

idx = nearest_idx3(outputIT(1:end-1), 0);
idx = nearest_idx3(outputIT(1:end-1), 0);

pre_stim_mean = mean(outputGau(1:idx)); 
post_stim_mean = mean(outputGau(idx:end)); 


if ~(max(mean_S_gau)) ==0
    ylim([0 max(mean_S_gau)])
    if size(t,2) > 1
        rectangle('position', [0 1 abs(mode(t(:,2)-t(:,1)))  max(mean_S_gau)], 'facecolor', [cfg.evt_color 0.5], 'edgecolor', [cfg.evt_color 0.5])
    else
        rectangle('position', [0 1 0.001  max(mean_S_gau)], 'facecolor', [cfg.evt_color 0.5], 'edgecolor', [cfg.evt_color 0.5])
    end
end

