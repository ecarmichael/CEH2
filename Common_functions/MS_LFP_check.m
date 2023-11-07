function MS_LFP_check(csc, offset, time_flag)
%% MS_LFP_check: plots all of the LFP channels with some offset.
%
%
%
%    Inputs:
%    - csc: [struct]     contains LFP signals in the TSD format.
%
%    - offset: [double]  how far to offset the LFP channels. default is
%    1000
%
%    - time_flag: [logical]   subtracts the starting time to bring it to
%    zero. default 1
%    Outputs:
%    -
%
%
%
%
% EC 2023-10-31   initial version
%
%
%
%% initialize

if nargin < 2
    offset = .001;
    time_flag = 1;
end
if nargin < 3
    time_flag = 1;
end

%%
figure(980)
clf
hold on
c_ord = linspecer(size(csc.data,1)+8);
c_ord(floor((size(csc.data,1)+4)/2):ceil((size(csc.data,1)+4)/2)+3,:) = [];

for ii = 1:size(csc.data,1)
    if time_flag
        plot(csc.tvec - csc.tvec(1), csc.data(ii,:)+offset*ii, 'color', c_ord(ii,:), 'linewidth', 1 )
        
    else
        plot(csc.tvec, csc.data(ii,:)+offset*ii, 'color', c_ord(ii,:), 'linewidth', 1 )
    end
    tick_val(ii) = csc.data(ii,1)+offset*ii;
    
end

set(gca, 'ytick', tick_val, 'YTickLabel', csc.label)
legend(csc.label, 'Location', 'northeast', 'Orientation', 'horizontal')
if time_flag
    xlim([csc.tvec(1), csc.tvec(end)] -csc.tvec(1))
else
    xlim([csc.tvec(1), csc.tvec(end)])
    
end
