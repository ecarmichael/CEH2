function [spk_count_all, cfg_out] = MS_event_hist(cfg_in,S, csc, varargin)
%% MS_event_hist: 
%
%
%    Inputs: 
%     - cfg [struct]   configuration see the defaults below.  
%
%     - S [struct] timestamp data.  To convert ms.Binary to S struct use MS_Binary2TS.m . 
%
%     - csc [struct] contains the TSD data for csc data.  At minimum has a
%     csc.tvec(nSamples) and csc.data(nChan, nSamples).  If empty will not
%     be used for plots. 
%
%     - IV_in [struct] contains timestamp data for the 
%
%     [ALTERNATIVE]
%     - tstart: [nEvents]  start times
%
%     - tend: [nEvents]   end times.  
%
%    Outputs: 
%     - spk_counts [nCell x nBin x nEvents]  histogram of spikes per bin
%     for each neuron for each event.  Can be used to get event-based
%     response hisogram. Alternatively can be used to get the response of a
%     cell across events.  
%
%
% EC 2020-04-29   initial version 
%
%
%%  initialize

if nargin == 4
    IV_in = varargin{1};
elseif nargin ==5
    IV_in = iv(varargin{1},varargin{2});
else
    error('Requires 4 or 5 inputs')
end

cfg_def = [];
cfg_def.remove_key = 'backspace'; % key that will flag events for removal.
cfg_def.t_win = [0 0];  % window around the events 
cfg_def.plot_chan = 'LFP';
cfg_def.ms_field = 'Binary'; % which ms field to use.   Can also be 'RawTraces'
cfg_def.bin_s = 0.01; 

cfg = ProcessConfig(cfg_def, cfg_in);
%% get the histogram of activity. 
   
centers = IVcenters(IV_in); 
IV_new = iv([centers+cfg.t_win(1) centers+cfg.t_win(end)]);
% IV_in.tstart = IV_in.tstart+cfg.t_win(1); 
% IV_in.tend = IV_in.tend+cfg.t_win(2); 

    % loop events
    for iEvt = length(IV_new.tstart):-1:1
        tbins = IV_new.tstart(iEvt):cfg.bin_s:IV_new.tend(iEvt); % make t bins
        tbins_center = tbins(1:end-1)+cfg.bin_s/2; % get the centers of the t bins
        
        for iC = length(S.t):-1:1
            spk_count = [];
            spk_count = histc(S.t{iC},tbins); % get spike counts for each bin
            spk_count_all(iC,:, iEvt) = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.
        end
        % sum for each event.  Helpful for plotting. 
        evt_activity(iEvt,:) = sum(spk_count_all(:,:,iEvt)); 
    end


%% plot each event if you like.  

if cfg.plot
    figure(112)
    ax(1) = subplot(5,1,1);
    if size(csc.data,2) >1
        cfg_plot.target = cfg.plot_chan;
    end
    
    PlotTSDfromIV(cfg_plot, IV_in, csc);
    xlim([csc.tvec(1) csc.tvec(end)])
    ax(2) = subplot(5,1,2:4);
    cfg_rast = [];
    cfg_rast.LineWidth  = 2;
    cfg_rast.openNewFig = 0;
    h = MultiRaster(cfg_rast,S);
%     set(gca, 'linewidth', 2)
    ylim([1 length(S.t)])
    vline(centers)
    linkaxes(ax, 'x')
    
    for iEvt = 1:length(IV_in.tstart)
        subplot(5,1,1);
        xlim([IV_new.tstart(iEvt), IV_new.tend(iEvt)])
        title(num2str(iEvt))
        subplot(5,1,5)
        bar(evt_activity(iEvt,:))
        ylim([0 10])
        
        
        drawnow
        pause
    end
    
    
end
    
end % of function





