function [IV_out, removed_evts] =  MS_check_IV(cfg_in,csc,ms_in,  varargin)
%% MS_check_IV: cycle through events and plot the csc along with the ms activity
%
%
%
%    Inputs: 
%     - cfg [struct]   configuration see the defaults below.  
%
%     - csc [struct] contains the TSD data for csc data.  At minimum has a
%     csc.tvec(nSamples) and csc.data(nChan, nSamples)
%
%     - ms_in [struct] ms miniscope structure.  used for plotting. 
%
%     - IV_in [struct] contains timestamp data for the 
%
%     [ALTERNATIVE]
%     - tstart: [nEvents]  start times
%
%     - tend: [nEvents]   end times.  
%
%    Outputs: 
%     - IV_out  [struct]   contains updated IV with flagged events removed.
%     
%
%
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
cfg_def.t_win = [-0.01 0.025];  % window around the events
cfg_def.plot_chan = 'LFP';
cfg_def.ms_field = 'Binary'; % which ms field to use.   Can also be 'RawTraces'

cfg = ProcessConfig(cfg_def, cfg_in);

%% plot and check. 
figure(112)
ax(1) = subplot(3,1,1);
if size(csc.data,2) >1
    cfg_plot.target = cfg.plot_chan;
end

PlotTSDfromIV(cfg_plot, IV_in, csc);
xlim([csc.tvec(1) csc.tvec(end)])


ax(2) = subplot(3,1,2:3);
traces = 1:size(ms_in.(cfg.ms_field),2);
for iT  = traces
    hold on
    this_B = ms_in.(cfg.ms_field)(:, iT);
    this_B(this_B == 0) = NaN; 
    plot(ms_in.time/1000,    (this_B*(0.01*length(traces)))+iT, 'linewidth', 2)
    this_B =[]; 
end
xlim([ms_in.time(1)/1000, ms_in.time(end)/1000]);
ylim([traces(1) traces(end)]);
linkaxes(ax, 'x')

remove_idx = zeros(1,length(traces));
for iEvt = 1:length(IV_in.tstart)
    xlim([IV_in.tstart(iEvt)-0.05, IV_in.tend(iEvt)+0.25])
    title(num2str(iEvt))
    drawnow
    waitforbuttonpress;
    key_hit = get(gcf, 'CurrentKey');
    if strcmp(key_hit, cfg.remove_key)
        remove_idx(iEvt) = 1;
    end
    
end

% update
IV_out = IV_in;
if ~isempty(remove_idx)
    IV_out.tstart(logical(remove_idx)) = [];
    IV_out.tend(logical(remove_idx)) = []; 
end

removed_evts = find(remove_idx); 

