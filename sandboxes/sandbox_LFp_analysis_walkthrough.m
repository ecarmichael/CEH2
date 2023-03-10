%% Sleep analysis tutorial

%% add code base and go to some data

mvdm_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
CEH2_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2'; % where the multisite repo can be found

data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\dSubiculum\inProcess\MD3\MD3_2022_07_20_D21'; 

cd(data_dir)

%% load some data

% load the NLX events
evts = LoadEvents([]);

blocks  = MS_get_waze_blocks(evts);

% load the NLX Spikes
cfg_s = [];
cfg_s.getTTnumbers = 0; % turn this off

S = LoadSpikes(cfg_s);


% load the position
cfg_p = [];
cfg_p.convFact = [6.95 6.95];
pos = MS_LoadPos([]);

% load the meta data
meta = MS_Load_meta;


% load the LFP and EMG
cfg_csc = [];
cfg_csc.fc = {meta.goodCSC, meta.EMG}; 
cfg_csc.desired_sampling_frequency = 2000; 
csc = MS_LoadCSC(cfg_csc);

emg = csc;
emg.data = emg.data(2,:); 
emg.label = {csc.label{2}}; 
emg.cfg.hdr = emg.cfg.hdr(1);

csc.data = csc.data(1,:);
csc.label = {csc.label{1}}; 
csc.cfg.hdr = csc.cfg.hdr(1); 

%%  plot the entire session

figure(101)
% plot(csc.tvec, csc.data(1,:))

plot(pos.data(1,:), pos.data(2,:), '.') 


%% restrict 

evts_w = restrict(evts, blocks.W_maze(1), blocks.W_maze(2)); 
S_w = restrict(S, blocks.W_maze(1), blocks.W_maze(2)); 
pos_w = restrict(pos, blocks.W_maze(1), blocks.W_maze(2)); 
csc_w = restrict(csc, blocks.W_maze(1), blocks.W_maze(2)); 


evts_of = restrict(evts, blocks.OF(1), blocks.OF(2)); 
S_of = restrict(S, blocks.OF(1), blocks.OF(2)); 
pos_of = restrict(pos, blocks.OF(1), blocks.OF(2)); 
csc_of = restrict(csc, blocks.OF(1), blocks.OF(2)); 


evts_post = restrict(evts, blocks.Post_sleep(1), blocks.Post_sleep(2)); 
S_post = restrict(S, blocks.W_maze(1), blocks.Post_sleep(2)); 
pos_post = restrict(pos, blocks.Post_sleep(1), blocks.Post_sleep(2)); 
csc_post = restrict(csc, blocks.Post_sleep(1), blocks.Post_sleep(2)); 
emg_post = restrict(emg, blocks.Post_sleep(1), blocks.Post_sleep(2)); 



%% test restriction

% figure(101)
% plot(csc_w.tvec, csc_w.data(1,:))

% plot(pos_w.data(1,:), pos_w.data(2,:), '.') 

%% Detect sleep


tic
    hypno = dSub_Sleep_screener(csc_post, emg_post); 
toc

%% SWR 

