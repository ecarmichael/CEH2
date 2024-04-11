function MS_get_binary_rate_maps(pos, data)


if isfield(pos, 'type')
    pos_type = 'tsd'; 
else
    pos_type = 'ms';
end


if isfield(data, 'type')
    data_type = 'tsd'; 
else
    data_type = 'ms';
end





% set up bins
SET_xmin = floor(min(pos.data(2,:))/100)*100 ; SET_ymin = floor(min(pos.data(1,:))/100)*100; % set up bins
SET_xmax = ceil(max(pos.data(2,:))/100)*100; SET_ymax = ceil(max(pos.data(1,:))/100)*100;

% SET_xmin = 100; 
% SET_xmax = 400;
% SET_ymin = 200;
% SET_ymax = 700;
SET_xBinSz =15; SET_yBinSz = 15;


 
x_edges = SET_xmin:SET_xBinSz:SET_xmax;
y_edges = SET_ymin:SET_yBinSz:SET_ymax;

% gaussian kernal
kernel = gausskernel([1 1],1); % 2d gaussian in bins

% compute occupancy encode
occ_hist = histcn(pos.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
occ_hist = conv2(occ_hist,kernel,'same');
 
no_occ_idx = find(occ_hist < 15); % NaN out bins never visited
occ_hist(no_occ_idx) = NaN;
 
occ_hist = occ_hist .* (1/30); % convert samples to seconds using video frame rate (30 Hz)

