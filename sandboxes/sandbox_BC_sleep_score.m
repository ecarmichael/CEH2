%% sandbox_BC_sleep_scoring 

data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Bryan_DropBox\CHRNA2_NOVEL_OBJECT\raw_data\NOPR\BC1602\BC1602_2023_06_30_D1_hab_failed';

emg_chan = 'CSC1.ncs';
lfp_chan = 'CSC3.ncs';


%% load some data

cd(data_dir)
cfg = [];
cfg.fc = {lfp_chan}; 

csc = MS_LoadCSC(cfg); 

% same for EMG
cfg_emg = [];
cfg_emg.fc = {emg_chan};
emg = MS_LoadCSC(cfg_emg); 

% filter the emg  % didn't actually do this here, but might be useful. 
% cfg_emg = [];
% cfg_emg.threshold = 0;
% cfg_emg.f = [15 30]; 
% % cfg_emg.type = 'cheby1';
% % cfg_emg.display_filter = 1
% emg = FilterLFP(cfg_emg,emg);



% load the tracking

% pos = LoadPos([]);  % commented out since you will be using DLC and not
% the nvt file. 

% load the events
evts = LoadEvents([]);


% restrict the data only to the sleep phase

start_OF = evts.t{find(contains(evts.label, 'Starting Recording'))}(1); 
start_sleep = evts.t{find(contains(evts.label, 'Starting Recording'))}(2); 

end_OF = evts.t{find(contains(evts.label, 'Stopping Recording'))}(1); 
end_sleep = evts.t{find(contains(evts.label, 'Stopping Recording'))}(2); 

% report the durations, to make sure there arn't any extra recordings
fprintf('<strong>OF duration: %.2f mins</strong>\n', (end_OF - start_OF)/60); 
fprintf('<strong>Sleep duration: %.2f mins</strong>\n', (end_sleep - start_sleep)/60); 


% restrict the data top just the sleep phase. 

csc_s = restrict(csc, start_sleep, end_sleep);
emg_s = restrict(emg, start_sleep, end_sleep);



%% plot a bit of data for quality check

figure(1)
clf

ax(1) = subplot(2,1,1);
plot(csc_s.tvec, csc_s.data);
legend('HC LFP')

ax(2) = subplot(2,1,2);
plot(emg_s.tvec, emg_s.data - .001, 'r') % offset it a bit
legend('emg')

linkaxes(ax, 'x'); % locks the x axes of the subplots so if you zoom on one the other zooms in. 

% fix the annoying wide limits on the x axis
xlim([csc_s.tvec(1) csc_s.tvec(end)])





%% sleep state

% specify known wake times or periods to ignore. 
wake_t = [6020 6722 6884 7437 7740 7870 7940 8140 8234 8460 9030 9170 11450 12740 15831 16898 17290 17418 17562 18133]; 
wake_idx = nearest_idx(wake_t, csc_s.tvec);
wake_idx = reshape(wake_idx,2, length(wake_idx)/2)'; 

%% score the sleep. 
[hypno, csc_out, emg_out] = dSub_Sleep_screener(csc_s, emg_s, wake_idx);  % can add in 'wake_idx' as the last input. 