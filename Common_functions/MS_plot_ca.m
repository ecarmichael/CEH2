function MS_plot_ca(cfg_in, ms_data_in, to_plot)
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
%     - to_plot: [string] shortcut to plot a specific data field. replaces
%     cfg.Ca_type; 
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

if nargin == 3
    cfg_in.Ca_type = to_plot; 
end


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
cfg_def.plot_type = '2d'; % '2d' or '3d'
cfg_def.offset = 1; 
cfg_def.rescale = []; 
cfg_def.width = 1.5; 
cfg_def.x_zoom = []; % where to zoom in on the x_axis for each plot.
cfg_def.view = [0 45]; 
cfg_def.colors = []; % can be used to pass a colour set. 
cfg_def.label = 'file_names';
cfg_def.saveas = []; % if empty don't save the images.  Can be '.fig' (matlab fig) or other known saveas format.  I like '.png'.

cfg = ProcessConfig(cfg_def, cfg_in);

if n_cells < 30
    cfg.Ca_chan = 1:n_cells;
end

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
            if strcmp(cfg.rescale, 'zscore')
                plot(time_in2*0.001, zscore(ms_data_in.(cfg.Ca_type)(:,cfg.Ca_chan(iC)))+iC*cfg.offset, 'color', c_ord(iC,:), 'linewidth', cfg.width)
                tick_val(iC) = median(zscore(ms_data_in.(cfg.Ca_type)(:,cfg.Ca_chan(iC)))+iC*cfg.offset);
            elseif strcmp(cfg.rescale, 'max')
                plot(time_in2*0.001, ms_data_in.(cfg.Ca_type)(:,cfg.Ca_chan(iC))./max(ms_data_in.(cfg.Ca_type)(:,cfg.Ca_chan(iC)))+iC*cfg.offset, 'color', c_ord(iC,:), 'linewidth', cfg.width)
                tick_val(iC) = median(ms_data_in.(cfg.Ca_type)(:,cfg.Ca_chan(iC))./max(ms_data_in.(cfg.Ca_type)(:,cfg.Ca_chan(iC)))+iC*cfg.offset);
            else
                plot(time_in2*0.001, ms_data_in.(cfg.Ca_type)(:,cfg.Ca_chan(iC))+iC*cfg.offset, 'color', c_ord(iC,:), 'linewidth', cfg.width)
                tick_val(iC) = median(ms_data_in.(cfg.Ca_type)(:,cfg.Ca_chan(iC))+iC*cfg.offset);
            end
            tick_label{iC} = cfg.Ca_chan(iC);
        end
        % 3d
    case '3d'
        for iC = 1:length(cfg.Ca_chan)
            if strcmp(cfg.rescale, 'zscore')
                plot3(time_in2*0.001, repmat(iC,size(ms_data_in.(cfg.Ca_type),1),1), zscore(ms_data_in.(cfg.Ca_type)(:,cfg.Ca_chan(iC))), 'color', c_ord(iC,:), 'linewidth', cfg.width)
            else
                plot3(time_in2*0.001, repmat(iC,size(ms_data_in.(cfg.Ca_type),1),1), ms_data_in.(cfg.Ca_type)(:,cfg.Ca_chan(iC)), 'color', c_ord(iC,:), 'linewidth', cfg.width)
            end
        end
        view(cfg.view)
end

xlim([time_in2(1)*0.001 time_in2(end)*0.001]);
y_lim = ylim;
ylim([0 y_lim(2)]); 

%         ax = [];
if ~isempty(cfg.x_zoom)
    
    xlim(cfg.x_zoom);
end

switch cfg.plot_type
    case '2d'
        [~, sort_idx] = sort(tick_val); % get past cases where the signal is lower than the offset previous. 
        set(gca, 'ytick', tick_val(sort_idx), 'yticklabel', tick_label(sort_idx));
        ylabel('cell id')
end
end % end function