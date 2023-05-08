%% sandbox_Trial_Wavelet_spec


%go to some data

data_dir = '/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Bryan_DropBox/ChRna2_sleep/Incoming/BC1743_sampling_light_pulses'

ft_dir = '/home/williamslab/Documents/Github/fieldtrip';

cd(data_dir);


%% load the data

cfg_csc = [];
cfg_csc.fc = {'CSC2.ncs'};

csc = MS_LoadCSC(cfg_csc);

% same but for EMG
cfg_emg= [];
cfg_emg.fc = {'CSC1.ncs'};

emg = MS_LoadCSC(cfg_emg);


%% load the events and find the light trigger times

evt = LoadEvents([]);


start_t_og = evt.t{(strcmp('Starting Recording', evt.label))}; 
start_t = evt.t{(strcmp('Starting Recording', evt.label))}; 

stop_t = evt.t{(strcmp('Stopping Recording', evt.label))}; 

pulse_t = evt.t{(strcmp('TTL Input on AcqSystem1_0 board 0 port 3 value (0x0020).', evt.label))}; % get the TTL times for the pulse. note that there is also an event for a TTL with the identifier (0x0000). This is the universal event logger.

% get the first pulse in each train
pdiff = diff(pulse_t);
pulse_start = pulse_t([0 find(pdiff > 5)]+1); % find all the pulses with a 5s gap before the previous pulse. Also include the first pulse as a 0 since we add 1 to correct for the diff. 

% find the short pulses
p_start_idx = [0 find(pdiff > 5)]+1; 

for ii = length(p_start_idx):-1:1
   p_type(ii)  = pulse_t(p_start_idx(ii)+1) - pulse_t(p_start_idx(ii)); 
end

pulse_fast = pulse_start(p_type < 1); % start of each pulse train with a pulse duration < 1; In the future you can you different cables and NLX TTL ports to have different event codes. This is just a hack for now. 
pulse_long = pulse_start(p_type >= 1); 


% correct time to start at 0
pulse_fast = pulse_fast - start_t_og;
pulse_long = pulse_long - start_t_og; 
pulse_start = pulse_start - start_t_og; 
pulse_t = pulse_t - start_t_og;
stop_t = stop_t - start_t_og;
start_t = start_t - start_t_og;


csc.tvec = csc.tvec - start_t_og;
emg.tvec = emg.tvec - start_t_og;

%% check the events make sense


figure(1010)
clf
plot(csc.tvec, csc.data(1,:))
hold on

xline(start_t,  '--k', 'start')
xline(stop_t, '--k', 'stop')

% for ii  = 1:length(pulse_t)
    vline(pulse_t, '-r');
        vline(pulse_start, '-b');

% end


%% get an event triggered average of the LFP 

win_s = .5; % window size in time
ETA = [];
for ii = length(pulse_start):-1:1
    
    ETA(ii,:) = csc.data(1,nearest_idx(pulse_start(ii) - win_s, csc.tvec):  nearest_idx(pulse_start(ii) + win_s, csc.tvec));
    
%     ETA(ii,:) = this_data.data(1,:);

end

% same thing but get random times as a shuffle comparison (takes some time)
Shuff_ts = csc.tvec(nearest_idx(csc.tvec(1)+win_s, csc.tvec): (nearest_idx(csc.tvec(end), csc.tvec))); % get a range of acceptibale times to pull from (cut off window size on each tend)

nShuff = 500;  % number of shuffles
ETA_shuff = [];

for iS = nShuff:-1:1
    this_shuff = randsample(Shuff_ts,1);
    ETA_shuff(iS,:) = csc.data(1,nearest_idx(this_shuff - win_s, csc.tvec):  nearest_idx(this_shuff + win_s, csc.tvec));
    
%     ETA(ii,:) = this_data.data(1,:);

end


%% plot the ETA
t_bins = -win_s:1/(csc.cfg.hdr{1}.SamplingFrequency):win_s;
figure(10)
clf
hold on
plot(t_bins, mean(ETA)); 
plot(t_bins, mean(ETA_shuff), '--', 'color', [.7 .7 .7]); 
xline(0, 'r')
xlabel('time from first pulse (s)')

legend({'ETA', 'Shuffle'})
SetFigure([], gcf); % makes figures nice for printing

%%  compute an event triggered spectrogram using FieldTrip

% uncomment these to add FT
% addpath(ft_dir) % add fieldtrip to the path. 
% ft_defaults % run this to set it up. 

figure
Triggered_Spec_FT(csc, pulse_long, [], 2:.5:120, [-5 -1], [-5 20])
xlim([-1 20])
vline([2:2:20])

figure
Triggered_Spec_FT(csc, pulse_fast, 'Short', 2:.5:120, [-5 -1], [-5 5])
xlim([-1 5])
vline([.1:.1:5])



