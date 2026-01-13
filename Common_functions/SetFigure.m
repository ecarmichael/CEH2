function h = SetFigure(cfg_in, h, sqr)
%  SetFigure will set the properties of figre "h" to the standards for EC
%  figure.
%          Inputs:
%           - cfg_in: [struct]  configuration parameters.  defaults:
%                 cfg_def.ft_size = 18;
%                 cfg_def.font = 'helvetica';
%                 cfg_def.grid = 'off';
%                 cfg_def.resize = 1;
%
%           - h : figure handle
%          Outputs:
%           - h : figure handle
% EC - 2016-10-05

%% defaults
if nargin < 3
    sqr = 0;
end

cfg_def.ft_size = 18;
cfg_def.font = 'helvetica';
% cfg_def.fontweight = 'normal';
cfg_def.grid = 'off';
cfg_def.resize = 1;
cfg_def.re_line = 0;

cfg = ProcessConfig2(cfg_def, cfg_in);


figure(h)
set(gcf,'windowstyle','normal');
set(h,'PaperPositionMode','auto')
set(gca,'DefaultTextFontSize',cfg.ft_size)
set(0, 'DefaulttextInterpreter', 'none')


if sqr == 1
    set(gcf,'units','normalized','outerposition',[.3 .1 .6 .8])
    
elseif cfg.resize == 1
    %     set(gcf, 'position', [600 50 560*1.4 420*1.4]);
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    %      set(gcf,'units','centimeters','outerposition',[0 0 10 10])
end
% H = get(gcf, 'children');
H = findobj(gcf,'type','axes','-not','Tag','legend','-not','Tag','Colorbar');

for iH = 1:length(H)
    set(H(iH), 'fontsize', cfg.ft_size, 'fontname', cfg.font, 'TickDir', 'out', 'fontweight', 'normal')
end

if cfg.re_line ~= 0
    for iH = 1:length(H)
        
        set(H(iH), 'linewidth', cfg.re_line)
    end
    
    
end

g = get(gca, 'XLabel');
set(g, 'fontsize', cfg.ft_size);
g = get(gca, 'YLabel');
set(g, 'fontsize', cfg.ft_size);
g = get(gca, 'title');
set(g, 'fontsize', cfg.ft_size);
box off



