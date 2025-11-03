function [d_ratio, p_ratio] =  MS_set_aspect_ratio(ax, d_scale, p_scale)
% rescales the axis ratios. ax is the axis, r_scale is the [x, y, z] scaling factors. A nice default is [.8 1 1];
if nargin < 1
    ax = gca;
    d_scale = [1 .67 1];
    p_scale = [1 .8  1];

elseif nargin < 2
    d_scale = [1 .67  1];
    p_scale = [1 .8  1];
elseif nargin < 3
    p_scale = [1 .8  1];
end

if isempty(ax)
    ax = gca;
end
y_lim = ylim;

set(ax, 'DataAspectRatio', get(gca, 'DataAspectRatio').*d_scale);
set(ax, 'PlotBoxAspectRatio', get(gca, 'PlotBoxAspectRatio').*p_scale);

ylim(y_lim)
d_ratio = get(ax, 'DataAspectRatio');
p_ratio = get(ax, 'PlotBoxAspectRatio');
