function MS_set_fig_pdf(h)
%% MS_set_fig_pdf:
%
%
%
%    Inputs: 
%    -
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2023-09-25   initial version 
%
%
%
%% initialize

if nargin < 1
    h = gcf; 
end

ft_size = 14; 

figure(h)
set(gcf,'windowstyle','normal');
set(h,'PaperPositionMode','auto')
set(gca,'DefaultTextFontSize',ft_size)
set(0, 'DefaulttextInterpreter', 'none')
% if cfg.resize == 1
    set(gcf, 'position', [600 50 560*1.4 420*1.4]);
% end
% H = get(gcf, 'children');
H = findobj(gcf,'type','axes','-not','Tag','legend','-not','Tag','Colorbar');

for iH = 1:length(H)
    set(H(iH), 'fontsize', ft_size, 'fontname', 'helvetica', 'TickDir', 'out', 'fontweight', 'normal')
end
box off
g = get(gca, 'XLabel');
set(g, 'fontsize',ft_size);
g = get(gca, 'YLabel');
set(g, 'fontsize', ft_size);
g = get(gca, 'title');
set(g, 'fontsize', ft_size);