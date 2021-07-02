%% dSub NLX screening sandbox


% this is a simple sandbox for loading, and checking the spatial tuning
% properties of cells using NLX


%% add some code to the path

 addpath(genpath('C:\Users\williamslab\Documents\github\vandermeerlab\code-matlab\shared'));

 addpath(genpath('C:\Users\williamslab\Documents\github\CEH2')); 

 data_dir = 'C:\Users\williamslab\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\dSubiculum\inProcess\M23_2021-07-02_OF';


cd(data_dir); % go to the data folder specified above
%% load dada
cfg.getTTnumbers = 0;
S = LoadSpikes(cfg);


pos = LoadPos([]);

% interpolate the spikes to match the time vector
spk_x = interp1(pos.tvec,pos.data(1,:),S.t{1},'linear');
spk_y = interp1(pos.tvec,pos.data(2,:),S.t{1},'linear');

%% plot the position

figure(101)

plot(pos.data(1,:), pos.data(2,:), '.', 'color', [0.8 0.8 0.8]);
hold on


S_idx = nearest_idx(pos.tvec, S.t{1}); 


plot(spk_x,spk_y, '.r')
axis off
