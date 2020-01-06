function Add_XY_scale(cfg_in)
%% Creates an appropriote X and Y scale line in the corner of a figure.  To run using defaults have a figure open and run "Add_XY_scale([])"
%
%
%  !!!!!CURRENTLY IN HACK MODE FOR FIXING FIGURES IN MULTISITE PAPER!!!!!
%
% Inputs:
%     - cfg_in [struct] contains adjustable aparmeters.  (TBD!!!)
%
%
%
%% Set defaults
cfg_def = [];
cfg_def.fontsize = 14;
cfg_def.hard_x = [];
cfg_def.hard_y_factor = 1;
cfg_def.hard_y =.5;
cfg_def.hard_y_factor = 1000; % LoadCSC is typically in 'V' so this is 1 'mV';
cfg_def.use_ticks = 1; % uses diff between ticks as the scale bar.
cfg_def.location = 'lower'; % where to put the scale bar. can be "upper" or "lower"
cfg_def.color = 'k'; % scale bar color;
cfg_def.linewidth = 2; % scale bar width;
% x defaults
cfg_def.xfactor = 1; % what fraction of the x axis range do you want. example Seconds to Ms use 1000.
cfg_def.xround = 1; % do you want to round the x scale? for example is data is in seconds and you want it in 100ms use:
cfg_def.xlabel = 'ms';
% y defaults
cfg_def.yfactor = 1000; % same for Y
cfg_def.yround = 1; % if true round to nearest whole
cfg_def.ylabel = 'mV';


cfg = ProcessConfig(cfg_def, cfg_in);


%% Get the X scale
xtick = get(gca, 'xtick');
if cfg.hard_x % use a hardcoded x value
    out_xtick = cfg.hard_x / cfg.hard_x_factor;
else
    if cfg.use_ticks
        scaled_xtick = ((xtick(end) - xtick(1))/(length(xtick)-1))*cfg.xfactor;
        fprintf('scaling X by factor %1.0f\n', cfg.xfactor)
        
    end
    if cfg.xround
        out_xtick = round(scaled_xtick,1);
    else
        
        out_xtick = scaled_xtick;
    end
end
%% Get the Y scale
ytick = get(gca, 'ytick');

if cfg.hard_y % use a hardcoded y value
    out_ytick = cfg.hard_y / cfg.hard_y_factor;
else
    
    
    if cfg.use_ticks
        scaled_ytick = ((ytick(end) - ytick(1))/length(ytick))*cfg.yfactor;
        fprintf('scaling Y by factor %1.0f\n', cfg.yfactor)
        
        if cfg.yround
            out_ytick = round(scaled_ytick,1);
        else
            out_ytick = scaled_ytick;
        end
    end
end
%% add to plot
fprintf('\n Xscale = %1.2f \n Yscale = %1.2f\n', out_xtick, out_ytick)

hold on
if strcmp(cfg.location, 'lower') || strcmp(cfg.location, 'Lower')
    y_offset = ytick(1) + (ytick(end) - ytick(1))/20;
elseif strcmp(cfg.location, 'upper') || strcmp(cfg.location, 'Upper')
    y_offset = ytick(end) - (ytick(end) - ytick(1))/20;
else
    error('unknown location. cfg.location must be "upper" or "lower"')
end


line([xtick(end-2) xtick(end-1)], [y_offset y_offset], 'color', cfg.color, 'linewidth', cfg.linewidth)

% same for ybar
line([xtick(end-1) xtick(end-1)], [y_offset y_offset+out_ytick], 'color', cfg.color, 'linewidth', cfg.linewidth)

if cfg.hard_x % use a hardcoded y value
    text(xtick(end-2)+((xtick(end-1)- xtick(end-2))/2), y_offset - (ytick(end) - ytick(1))/20, [num2str(cfg.hard_x), cfg.xlabel], 'fontsize', cfg.fontsize, 'fontname', 'helvetica')
    
else
    text(xtick(end-2)+((xtick(end-1)- xtick(end-2))/2), y_offset - (ytick(end) - ytick(1))/20, [num2str(out_xtick), cfg.xlabel], 'fontsize', cfg.fontsize, 'fontname', 'helvetica')
end


if cfg.hard_y % use a hardcoded y value
    
    text(xtick(end-1)+(xtick(end) - xtick(1))/50, median([y_offset y_offset+out_ytick]), [num2str(cfg.hard_y), cfg.ylabel], 'fontsize', cfg.fontsize, 'fontname', 'helvetica')
else
    text(xtick(end-1)+(xtick(end) - xtick(1))/50, median([y_offset y_offset+out_ytick]), [num2str(out_ytick), cfg.ylabel], 'fontsize', cfg.fontsize, 'fontname', 'helvetica')
end
% median([y_offset y_offset+out_ytick]), '

% reset the axis
% ylim(


