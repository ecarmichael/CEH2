function MS_plot_ca(cfg_in, ms_data_in)
%% MS_plot_ca_nlx: makes a simple plot of calicum traces. Can be either segmented (recording blocks) or one continous signal.
%
%
%
%    Inputs:
%     - cfg [struct] contains configuration paramters [see below]
%     - ms_data_in
%           if continuous, this is the ms.RawTraces or ms.FiltTraces [time x cell]
%           if segmented, this is a cell array with each containing
%           RawTraces or Filttraces [nSeg cells]
%
%
%    Outputs:
%     - h : plot handle.
%
%       ToDO:
%           - switch based on ms_in being cell or not
%           - use csc channel labels to pick a channel when more than one
%           exist.
%
%
% EC 2020-01-20   initial version
%
%
%% initialize

% if the ms_data in is continous or segmented (in cells) set the type.
if ~iscell(ms_data_in.RawTraces)
    type = 'continuous';
    n_cells = size(ms_data_in.RawTraces,2);
    n_seg = 1;
else
    type = 'segmented';
    n_cells = size(ms_data_in.RawTraces{1},2);
    n_seg = length(ms_data_in.RawTraces); % how many recording segments.
end

%% set up defaults.
cfg_def =[];
cfg_def.Ca_type = 'RawTraces'; % can be either 'RawTraces' or FiltTraces' (maybe others?)
cfg_def.Ca_chan = 1:floor(n_cells/20):n_cells; % get a subset of cells.
cfg_def.plot_type = '3d'; % '2d' or '3d'
cfg_def.x_zoom = []; % where to zoom in on the x_axis for each plot.
cfg_def.view = [0 45]; 
cfg_def.colors = []; % can be used to pass a colour set. 
cfg_def.label = 'file_names';
cfg_def.saveas = []; % if empty don't save the images.  Can be '.fig' (matlab fig) or other known saveas format.  I like '.png'.

cfg = ProcessConfig(cfg_def, cfg_in);

%% make the plots
if isempty(cfg.colors)
    c_basic = linspecer(10); % basic color range
    c_basic(5:6,:) = [];
    c_ord = repmat(c_basic, (ceil(size(cfg.Ca_chan,2)/length(c_basic))),1); % repeat the colors range for better visibility;
    % c_ord = linspecer(length(cfg.Ca_chan)); % nice colours.
else
    c_ord = cfg.colors;
end

hold on
time_in2 = ms_data_in.time - ms_data_in.time(1);
switch cfg.plot_type
    % 2d
    case '2d'
        for iC = 1:length(cfg.Ca_chan)
            
            plot(time_in2*0.001, ms_data_in.(cfg.Ca_type)(:,iC)+iC*0.5, 'color', c_ord(iC,:))
            tick_val(iC) = median(ms_data_in.(cfg.Ca_type)(:,iC)+iC*0.5); 
            tick_label{iC} = num2str(iC); 
        end
        % 3d
    case '3d'
        for iC = 1:length(cfg.Ca_chan)
            plot3(time_in2*0.001, repmat(iC,size(ms_data_in.(cfg.Ca_type),1),1), ms_data_in.(cfg.Ca_type)(:,iC), 'color', c_ord(iC,:))
        end
        view(cfg.view)
end

xlim([time_in2(1)*0.001 time_in2(end)*0.001]);
%         ax = [];
if ~isempty(cfg.x_zoom)
    
    xlim(cfg.x_zoom);
end

switch cfg.plot_type
    case '2d'
        set(gca, 'ytick', tick_val, 'yticklabel', tick_label);
        ylabel('cell id')
end
end % end function