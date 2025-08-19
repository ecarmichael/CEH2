function [dist, cdist] = MS_path_dist_spd(x, y, t, plt_flag, xlims, ylims)
%% MS_path_dist: calculates the cumulative distance along a path using the trapz method
%
%
%
%   Inputs:
%       x: [array]  x positions
%
%       y: [array]  t positions
%
%       t: [array]  time
%
%    Optional:       
%       plt_flag: [logical] do you want a plot of the path and cumulative dist? 
%       xlims: [2 x 1]  limits for inclusion (x)
%
%       ylims: [2 x 1]  limits for inclusion (y)


%% init

error('DO NOT USE> overestimates distance due to jumps in speed. ')

if nargin < 4
    plt_flag = 0; 
    xlims = [min(x) max(x)];
    ylims = [min(y) max(y)];
elseif nargin <5
    xlims = [min(x) max(x)];
    ylims = [min(y) max(y)];
end




%% remove outliers and interpolate

o_idx = x < xlims(1) | x > xlims(2); 
x_i = x; 
x_i(o_idx) =nan; 

o_idx = y < ylims(1) | y > ylims(2); 
y_i = y; 
y_i(o_idx) =nan; 

p_i(1,:) = fillmissing(x_i, 'nearest');
p_i(2,:) = fillmissing(y_i, 'nearest');

%% get the distance and cumulative distance vector
% based on: https://www.mathworks.com/matlabcentral/answers/2096886-what-is-the-best-way-of-calculating-the-path-length-of-a-freely-moving-tracked-animal


x = p_i(1,:) - min(p_i(1,:)); 
y = p_i(2,:) - min(p_i(2,:)); 

dxdt = gradient(x, t); 
dydt = gradient(y, t);

dist = trapz(t, sqrt(dydt.^2 + dxdt.^2));

cdist = cumtrapz(t, sqrt(dydt.^2 + dxdt.^2));









