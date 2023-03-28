fs = 32000; 
chanMap = 1:64;
chanMap0ind = chanMap - 1;
connected = true(64, 1);
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

% kcoords = reshape(repmat(1:4, 1,16), 1,64);

save('Shuttle_16TTChanMap.mat', ...
    'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs')