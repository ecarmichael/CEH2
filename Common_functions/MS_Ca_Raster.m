function h = MS_Ca_Raster(data_in, tvec, linesize, color_in)
%% MS_Ca_Raster:
%
%
%
%    Inputs: 
%    - data_in: [nChan x nSample] matrix with binary events.  
%
%    - tvec: [1 x nSample] time vector [optional]
%
%
%
%    Outputs: 
%    - h [handle] 
%
%
%
%
% EC 2021-03-27   initial version 
%
%
%
%% initialize

if nargin < 2 || isempty(tvec)
    tvec = 1:size(data_in,2);
    fprintf('No tvec give, using data length (%i samples) as vector\n', size(data_in, 2))
    linesize = 2; 
    color_in = []; 
elseif nargin <3
    linesize = 2; 
    color_in = [];
elseif nargin < 4
    color_in = []; 
end

if size(tvec,2) ~= size(data_in, 2)
    error('data_in and tvec differ in size')
end
%%
nChan = size(data_in, 1); 
% set colours
if isempty(color_in)
c_ord = linspecer(nChan+1); % add one more for background
else
    c_ord = repmat(color_in, nChan+1,1);
end
j_ord = jet(100); % set the set 

hold on
set(gca, 'YDir','reverse'); % flip y axis order
xlim([tvec(1) tvec(end)]);
ylim([1, nChan])

for iChan = nChan:-1:1

    % convert data_in into indicies. this is essential to it running without crashing everything.
    idx = find(data_in(iChan,:) ~= 0);
    
    % plot the line using only the indices. 
    if ~isempty(idx)
        line([tvec(idx); tvec(idx)], [iChan-2 iChan+1], 'color', c_ord(iChan+1,:), 'linewidth', linesize)
    end
end

set(gca, 'color', j_ord(1,:)); %set background color. 
