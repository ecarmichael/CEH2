%% SWR calcium sandbox

addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));

ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1252\2021_12_16_pv1252_MZD3'; 

lfp_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\LFP\2021_12_16_pv1252_MZD3_LFP';

%% find SWR for one session
cd(lfp_dir)
cfg = []; 
cfg.fc = {'CSC6.ncs', 'CSC7.ncs'};
cfg.desired_sampling_frequency = 2000;
csc = MS_LoadCSC(cfg);

pos = LoadPos([]);

%% restrict to track periods
evt = LoadEvents([]);

rec_s_idx  = find(contains(evt.label, 'Starting Recording')); 
rec_e_idx  = find(contains(evt.label, 'Stopping Recording')); 

for ii= length(evt.t{rec_s_idx}):-1:1
    rec_len(ii) = evt.t{rec_e_idx}(ii) - evt.t{rec_s_idx}(ii);
    fprintf('Rec # %.0f: %2.2fmins\n', ii, rec_len(ii)/60)
end

[sort_dur, sort_idx] = sort(rec_len); 

keep_rec = sort_idx(end-2:end); 
keep_rec = sort(keep_rec); 
pre_sleep = [evt.t{rec_s_idx}(keep_rec(1)), evt.t{rec_s_idx}(keep_rec(1))+60*60*2]; 
maze = [evt.t{rec_s_idx}(keep_rec(2)) evt.t{rec_e_idx}(keep_rec(2))]; 
post_sleep = [evt.t{rec_s_idx}(keep_rec(3)), evt.t{rec_s_idx}(keep_rec(3))+60*60*2]; 


pre_csc = restrict(csc, pre_sleep(1), pre_sleep(2)); disp((pre_csc.tvec(end) - pre_csc.tvec(1))/60)
maze_csc = restrict(csc, maze(1), maze(2)); disp((maze_csc.tvec(end) - maze_csc.tvec(1))/60)
maze_pos = restrict(pos,  maze(1), maze(2)); maze_pos.data = (maze_pos.data/.666)/10; maze_pos.units = 'cm'; 
post_csc = restrict(csc, post_sleep(1), post_sleep(2));disp((post_csc.tvec(end) - post_csc.tvec(1))/60)

nan_idx = []; 

this_csc = maze_csc; 
this_csc_swr = this_csc;
this_csc_swr.data = maze_csc.data(1,:); 
this_csc_swr.label = maze_csc.label(1); 
this_csc_swr.cfg.hdr = []; 
this_csc_swr.cfg.hdr{1} = maze_csc.cfg.hdr{1};

this_csc.data = maze_csc.data(2,:); 
this_csc.label = maze_csc.label(2); 
this_csc.cfg.hdr = [];
this_csc.cfg.hdr{1} = maze_csc.cfg.hdr{2};
%% select indicies for exclusion

maze_linspd = getLinSpd([],maze_pos); 
maze_linspd.data = smooth(maze_linspd.data, ceil(1/mode(diff(maze_pos.tvec))))'; 

speed_int = interp1(maze_linspd.tvec, maze_linspd.data, maze_csc.tvec);

move_idx = speed_int > 10; 
sat_idx = abs(maze_csc.data) == max(abs(maze_csc.data)); % remove saturations. if non this will just remove the max point. 

nan_idx = [sat_idx | move_idx' ]; 

%% get the contrast filters and their ratios

cfg_con_f = [];
cfg_con_f.threshold = 0;
cfg_con_f.f = [6 12]; cfg_con_f.type = 'fdesign';
%     cfg_con_f.display_filter = 1
csc_th = FilterLFP(cfg_con_f,this_csc);
% csc_th.data = smooth(csc_th.data, csc.cfg.hdr{1}.SamplingFrequency*2);
csc_th.data = smooth(abs(hilbert(csc_th.data)), this_csc.cfg.hdr{1}.SamplingFrequency);

cfg_con_f = [];
cfg_con_f.threshold = 0;
cfg_con_f.f = [1 4]; cfg_con_f.type = 'fdesign';
%     cfg_con_f.display_filter = 1
csc_delta = FilterLFP(cfg_con_f,this_csc);
% csc_delta.data = smooth(csc_delta.data, csc.cfg.hdr{1}.SamplingFrequency*2);
csc_delta.data = smooth(abs(hilbert(csc_delta.data)), this_csc.cfg.hdr{1}.SamplingFrequency);


ratio = this_csc;
ratio.data = ((csc_th.data ./ csc_delta.data)');


% z_ratio = zscore(csc_th.data./csc_delta.data);

if isempty(nan_idx)
nan_idx = (ratio.data >.5); 
else
%     nan_idx = interp1(ratio.tvec, ratio.
    nan_idx = nan_idx | (ratio.data >.5);
end

%% check the movement and Theta/Delta exclusion
figure(1001)
clf
subplot(2,1,1)
hold on
plot(maze_csc);
plot(maze_csc.tvec, nan_idx/1000)
xlim([maze_csc.tvec(1) maze_csc.tvec(end)])

subplot(2,1,2)
hold on
plot(maze_csc.tvec, nan_idx-1)
plot(this_csc.tvec, speed_int)
xlim([maze_csc.tvec(1) maze_csc.tvec(end)])

%% filter the LFP into the ripple band.
cfg_swr = [];
cfg_swr.check = 0; % plot checks.
cfg_swr.filt.type = 'butter'; %Cheby1 is sharper than butter
cfg_swr.filt.f  = [125 200]; % broad, could use 150-200?
cfg_swr.filt.order = 4; %type filter order (fine for this f range)
cfg_swr.filt.display_filter = 0; % use this to see the fvtool

% smoothing
cfg_swr.kernel.samples = csc.cfg.hdr{1}.SamplingFrequency/100;
cfg_swr.kernel.sd = csc.cfg.hdr{1}.SamplingFrequency/100;

% detection
% cfg_swr.artif_det.method = 'zscore';
% cfg_swr.artif_det.threshold = 8;
% cfg_swr.artif_det.dcn = '>';
% cfg_swr.artif_det.rm_len = .2;
cfg_swr.threshold = 1.5;% in sd
cfg_swr.method = 'zscore';
cfg_swr.min_len = 0.04; % mouse SWR: 40ms from Vandecasteele et al. 2014
cfg_swr.merge_thr = 0.02; %merge events that are within 20ms of each other.

% 
% cfg_swr.nan_idx = nan_idx; % where are any nans, say from excluding artifacts, other events...

% restrictions
cfg_swr.max_len = [];
cfg_swr.max_len.operation = '<';
cfg_swr.max_len.threshold = .1;

%                 cfg_swr.min_len = [];
%                 cfg_swr.min_len.operation = '<';
%                 cfg_swr.min_len.threshold = .2;
cfg_swr.nCycles = 5; % number of cycles
cfg_swr.nCycles_operation = '>='; % number of cycles

% variaence
% cfg_swr.var = [];
% cfg_swr.var.operation = '<';
% cfg_swr.var.threshold = 1;

SWR_evts = MS_get_LFP_events_sandbox(cfg_swr, this_csc);

%
%         cfg_max_len = [];
%         cfg_max_len.operation = '>';
%         cfg_max_len.threshold = 5;
%          SWR_evts = SelectIV(cfg_max_len,SWR_evts,'nCycles');
subplot(2,1,1)
% check quality. 
cfg_plot.display = 'tsd'; %'iv';
cfg_plot.title = 'var_raw';
PlotTSDfromIV(cfg_plot, SWR_evts, this_csc)

