fs = 30000; 

chanMap = [46 47 44 42 45 43 38 40 41 36 34 32 39 30 28 26 37 35 24 22 33 31 20 18 29 27 25 16 23 21 19 17, 14 15 12 10 13 11 6 8 9 4 2 0 7 62 60 58 5 3 57 54 1 63 52 50 61 59 56 48 55 53 51 49]';
chanMap0ind = chanMap - 1;
connected = true(64, 1);
% bad_chan = [3 26 28 31 32 33 34 37 43 44 45 46 53 54 63];
% connected(bad_chan) = 0; 
% xcoords =  reshape(repmat(5:5:20, 1,16), 1,64); 
% ycoords = [20 20 20 20 40 40 40 40 60 60 60 60 80 80 80 80 100 100 100 100 120 120 120 120 140 140 140 140 160 160 160 160, [ 20 20 20 20 40 40 40 40 60 60 60 60 80 80 80 80 100 100 100 100 120 120 120 120 140 140 140 140 160 160 160 160]+160]';
x = [20     0    40     20 ];%   25    20    30    25  45    40    50    45 65 60 70 65];
y = [0     20     20    40];%   20    25    25    30  40    45    45    50 60 65 65 70];
xcoords = [];
ycoords = [];
kcoords = []; 
for ii = 1:16
    kcoords = [kcoords, ii ii ii ii];
    xcoords = [xcoords, x+(500)*ii];
    ycoords = [ycoords, y+(500)*ii]; 
end
kcoords = kcoords';
xcoords = xcoords';
ycoords = ycoords';
% chanMap(bad_chan) = [];
% chanMap0ind(bad_chan) = [];
% connected(bad_chan) = [];
% xcoords(bad_chan) = [];
% ycoords(bad_chan) = [];
% kcoords(bad_chan) = [];
% kcoords = reshape(repmat(1:4, 1,16), 1,64);

save('Shuttle_16TT_intanChanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')

%% reshape the coords to match the channel labels. 
chan_idx = chanMap;
chanMap = [1:64]';
chanMap0ind = chanMap - 1;
connected = true(64, 1);
 bad_chan = [3 26 28 31 32 33 34 37 43 44 45 46 53 54 63];
connected(bad_chan) = 0; 
% xcoords_og = xcoords;


xcoords = xcoords(chan_idx+1); 
ycoords = ycoords(chan_idx+1); 
kcoords = kcoords(chan_idx+1); 

save('Shuttle_16TT_test.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')
