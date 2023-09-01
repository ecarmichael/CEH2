function MS_rescale_axis(x, y, scale)
%% MS_rescale_axis: adjusts the xlim and ylim by some percentage of the range. 
%
%
%
%    Inputs: 
%    - x: [1 x N]     x data
%    - y: [1 x N]     y data
%
%    - scale: [double] or [1 x 2   percentage to scale by some
%    percentage [x y]


%    example scale = [0.1] will
%    make everything 10% wider in all axes.   
%
%    Example 2: scale  = [.1 .3] will scale x by 10% of the x range and 30%
%    of the y range
%
%    Outputs: 
%    - none
%
%
%
%
% EC 2023-08-31   initial version 
%
%% initialize

if length(scale) == 1
    scale = [scale, scale]; 
end
%
%% 

xlim([min(x) - ((max(x) - min(x))*scale(1)), max(x) + ((max(x) - min(x))*scale(1))]);
ylim([min(y) - ((max(y) - min(y))*scale(2)), max(y) + ((max(y) - min(y))*scale(2))]);



