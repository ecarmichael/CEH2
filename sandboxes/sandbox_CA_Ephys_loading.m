%% CA2 + ephys sandbox

%% add paths

close all
restoredefaultpath
global PARAMS
os = computer; 

if ismac
    PARAMS.data_dir = '/Users/jericcarmichael/Documents/Williams_Lab/7_12_2019_PV1069_LTD5'; % where to find the raw data
    PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/Users/jericcarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found

elseif strcmp(os, 'GLNXA64')

    PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/7_12_2019_PV1069_LTD5'; % where to find the raw data
    PARAMS.inter_dir = '/home/ecarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/home/ecarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found

else
 disp('on a PC')
end


rng(11,'twister') % for reproducibility


% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_CEH2_dir));
cd(PARAMS.data_dir) % move to the data folder

% try the newer NLX loaders for UNIX
[~, d] = version;
if str2double(d(end-3:end)) >2014
    rmpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared/io/neuralynx')
    addpath(genpath('/Users/jericcarmichael/Documents/NLX_loaders_UNIX_2015'))
    disp('Version is greater than 2014b on UNIX so use updated loaders found here:')
    which Nlx2MatCSC
end

clear d os
%%  Laod some stuff
tic
load('ms.mat');
% load('SFP.mat');
TS = MS_Load_TS('timestamp.dat');

% [MS_ts.camNum,MS_ts.frameNum,MS_ts.sysClock,MS_ts.buffer1] = MS_Load_TS('timestamp.dat');

% %% make a video
% 
% for iframe =  1:size(SFP, 3)
%     
%     imagesc(SFP(:,:,iframe))
%     M(iframe) = getframe;
% end

toc
%% load MS timestamps
TS= MS_Load_TS('timestamp.dat');
TS2 = MS_Load_TS('timestamp2.dat');

% [TS, CAMNUM,FRAMENUM,SYSCLOCK,BUFFER1] = MS_Load_TS('timestamp.dat');

% if length(unique(CAMNUM)) 
% this_cam = 0;
%     for ii = unique(CAMNUM)'
% %         disp(ii)
%         this_cam = this_cam+1;
%         idx = CAMNUM == ii;
%     TS.CAMNUM{this_cam} = CAMNUM(idx);
%     TS.FRAMENUM{this_cam} = FRAMENUM(idx);
%     TS.SYSCLOCK{this_cam} = SYSCLOCK(idx);
%     TS.BUFFER1{this_cam} = BUFFER1(idx);
% 
%     end
   

%% GET nlx data
tic
cfg = [];
cfg.fc = {'CSC4.ncs'};%, 'CSC8.ncs'};
csc = LoadCSC(cfg);

cfg = [];
evt = LoadEvents(cfg);
evt.t{length(evt.t)+1} = unique(sort([evt.t{3} evt.t{4}]));
evt.label{length(evt.label)+1} = 'all_evt';

% get the recording start time from NLX header
if isfield(csc.cfg.hdr{1}, 'TimeCreated')
    NLX_start = csc.cfg.hdr{1}.TimeCreated; % find the creation time as a string
    if contains(NLX_start, ':')
        NLX_start = duration(str2double(strsplit(NLX_start(end-8:end),':'))); % pull out hours:mins:sec and convert to a time
    else
        NLX_start = duration([NLX_start(end-5:end-4) ':' NLX_start(end-3:end-2) ':' NLX_start(end-1:end)]); % pull out hours:mins:sec and convert to a time
    end
end

% identify major jumps in evts

 all_jumps = diff(evt.t{5}) > (mean(diff(evt.t{5}) +2*std(diff(evt.t{5}))));
 all_jumps(1) = 0; % correct for first jump;
 jump_idx = find(all_jumps ==1);
 rec_evt = [];
 if sum(all_jumps) > 0 && sum(all_jumps) <2
     fprintf('Jump found at time: %.0f\n', evt.t{5}(jump_idx))
     
     rec_evt{1} = restrict(evt, evt.t{5}(1), evt.t{5}(jump_idx)); % add one index to compensate for the diff. 
     rec1_csc = restrict(csc, evt.t{5}(1), evt.t{5}(jump_idx));
     
     rec_evt{2} = restrict(evt, evt.t{5}(jump_idx+1), evt.t{5}(end));
     rec2_csc = restrict(csc, evt.t{5}(jump_idx+1), evt.t{5}(end));
     
 elseif sum(all_jumps) >2
     
     for iJ = length(jump_idx):-1:1
         if iJ ==1
             rec_evt{iJ} = restrict(evt, evt.t{5}(1), evt.t{5}(jump_idx(iJ)));
%              rec_csc{iJ} = restrict(csc, evt.t{5}(1), evt.t{5}(jump_idx(iJ)));
         else
             rec_evt{iJ} = restrict(evt, evt.t{5}(jump_idx(iJ-1)), evt.t{5}(jump_idx(iJ))); 
%              rec_csc{iJ} = restrict(csc, evt.t{5}(jump_idx(iJ-1)), evt.t{5}(jump_idx(iJ))); 

         end
     end
 end
toc

%% if the TSs align with the evt then add it in as a subfield

for iE = 1:length(rec_evt)
    this_cam = 2;
    if length(rec_evt{iE}.t{5})== length(TS{iE}.system_clock{this_cam})
        disp(['TS and evts align! for Timestamp: ' TS{iE}.filename])
        TS{iE}.NLX_ts = rec_evt{iE}.t{5}';
        TS{iE}
    else
        warning(['TS and evt do not align for Timestamp: ' TS{iE}.filename])
    end
end

%% make a interpolated signal to see where things are missing
ts_norm = [0 ;TS{1}.system_clock{2}(2:end)]';
evt_norm = rec_evt{1}.t{5} - rec_evt{1}.t{5}(1); 


%% check length of TSs
disp('TS1')
for this_cam = 1:length(TS.system_clock)
fprintf('Number of Scope TS id: %.0f  =   %.0f  at %0.2f Hz\n',this_cam, length(TS.system_clock{this_cam}), 1/(median(diff(TS.system_clock{this_cam}(2:end)))*0.001))
end

    
disp('Rec1')
for this_evt = 3:length(rec1_evt.label) % correct for start and stop recording. 
    fprintf('Number of evt evts id: %.0f  =   %.0f at %0.2f Hz\n',this_evt, length(rec1_evt.t{this_evt}),1/(median(diff(rec1_evt.t{this_evt}))))
end

disp('TS2')
for this_cam = 1:length(TS2.system_clock)
fprintf('Number of Scope TS id: %.0f  =   %.0f  at %0.2f Hz\n',this_cam, length(TS2.system_clock{this_cam}), 1/(median(diff(TS2.system_clock{this_cam}(2:end)))*0.001))
end

disp('Rec2')
for this_evt = 3:length(rec2_evt.label) % correct for start and stop recording. 
    fprintf('Number of evt evts id: %.0f  =   %.0f at %0.2f Hz\n',this_evt, length(rec2_evt.t{this_evt}),1/(median(diff(rec2_evt.t{this_evt}))))
end


disp('All evt')
for this_evt = 3:length(evt.label) % correct for start and stop recording. 
    fprintf('Number of evt evts id: %.0f  =   %.0f at %0.2f Hz\n',this_evt, length(evt.t{this_evt}),1/(median(diff(evt.t{this_evt}))))
end

% fprintf('Number of NLX events: %.0f, Number of Scope TS: %.0f, Difference: %.0f\n', length(evt.t{this_evt}), length(TS.SYSCLOCK{this_cam}(2:end)), length(evt.t{this_evt}) -  length(TS.SYSCLOCK{this_cam}(2:end)));


%% give new TS to MS data 
% if things checkout
TS2.NLX_tvec{2} = rec2_evt.t{5}';
TS2.NLX_tvec{1} = interp(rec2_evt.t{5}, 2)';

%% for JISU DATA

peak_threshold =  (mean(diff(evt.t{5}) +2*std(diff(evt.t{5}))));
min_dist = 10;
[~, Rec_idx] = findpeaks(diff(evt.t{5}), 'minpeakheight',peak_threshold, 'minpeakdistance', min_dist);
fprintf(['\nDetected %.0f trigger transitions treating this as %.0f distinct recordings\n'], length(Rec_idx), length(Rec_idx)/2)


t_start = Rec_idx(1:2:end-1);
t_end = Rec_idx(2:2:end);


 rec_evt = []; 

for iE = 1:length(t_start)
    rec_evt{iE} = restrict(evt, evt.t{5}(t_start(iE)), evt.t{5}(t_end(iE)));
end
    
figure(1)
hold on
plot(diff(evt.t{5}), 'k')
hline(peak_threshold, '--r')
plot(Rec_idx, 100, '*k')

% break them into recording sessions
t_start = Rec_idx(1:2:end-1);
t_end = Rec_idx(2:2:end);


plot([t_start ; t_end]', [50 50], '-b')


%      for iJ = length(Rec_idx):-1:1
%          if iJ ==1
%              rec_evt{iJ} = restrict(evt, evt.t{5}(1), evt.t{5}(Rec_idx(iJ)));
% 
%          else
%              rec_evt{iJ} = restrict(evt, evt.t{5}(Rec_idx(iJ-1)), evt.t{5}(Rec_idx(iJ)));
%          end
%      end


%% check length of TSs

for iRec = 1:length(rec_evt)
    disp(['Rec ' num2str(iRec)])
    for this_evt = 3:length(rec_evt{iRec}.label) % correct for start and stop recording.
        fprintf('Number ofevts id: %.0f  =   %.0f samples at %0.2f Hz. Start at: %s for ~%.1fs \n',this_evt, length(rec_evt{iRec}.t{this_evt}),1/(median(diff(rec_evt{iRec}.t{this_evt}))),...
            char(NLX_start + minutes((rec_evt{iRec}.t{this_evt}(1)  - csc.tvec(1))/60)), length(rec_evt{iRec}.t{this_evt})/(1/(median(diff(rec_evt{iRec}.t{this_evt})))))
    end
end

for iRec = 1:length(TS)
    disp(['TS ' num2str(iRec)])
    fprintf('Number of Scope TS id: %.0f  =   %.0f  at %0.2fHz for %.f sec\n',iRec, length(TS{iRec}.system_clock{1}), 1/(median(diff(TS{iRec}.system_clock{1}(2:end)))*0.001),...
        length(TS{iRec}.system_clock{1})/ (1/(median(diff(TS{iRec}.system_clock{1}(2:end)))*0.001)))
    
end


disp('All evt')
for this_evt = 3:length(evt.label) % correct for start and stop recording. 
    fprintf('Number of evt evts id: %.0f  =   %.0f at %0.2f Hz\n',this_evt, length(evt.t{this_evt}),1/(median(diff(evt.t{this_evt}))))
end

% fprintf('Number of NLX events: %.0f, Number of Scope TS: %.0f, Difference: %.0f\n', length(evt.t{this_evt}), length(TS.SYSCLOCK{this_cam}(2:end)), length(evt.t{this_evt}) -  length(TS.SYSCLOCK{this_cam}(2:end)));


%% restrict data to first recording of the session
csc_r = restrict(csc, evt.t{1}(1), evt.t{2}(1));
evt_r = restrict(evt, evt.t{1}(1), evt.t{2}(1));


% if evt_r.t{3}(1) < evt_r.t{4}(1)
%     all_evts_r = [evt_r.t{3} evt_r.t{4}];
% elseif evt_r.t{3}(1) > evt_r.t{4}(1)
%     all_evts_r = [evt_r.t{4} evt_r.t{3}];
% elseif evt_r.t{3}(1) == evt_r.t{4}(1)
%     warning(['Event times for ' evt_r.label{3} ' are somehow equal to ' evt_r.label{4} '.  Check into this...'])
% end

%% correct for recording time (just to make things easier)
% for ii = 1:length(evt_r.t)
%     evt_r.t{ii} = evt_r.t{ii} - csc_r.tvec(1);
% end

    all_evts_r = unique(sort(evt.t{5}));

    % set a max cutoff_just for plotting
    

% all_evts_r = all_evts_r - csc_r.tvec(1);

% csc_r.tvec = csc_r.tvec - csc_r.tvec(1);

%% get some recording periods
peak_threshold = 5;
[~, Rec_idx] = findpeaks(diff(all_evts_r), 'minpeakheight',peak_threshold);
fprintf(['\nDetected %.0f trigger transitions treating this as %.0f distinct recordings\n'], length(Rec_idx), length(Rec_idx)/2)

% plot the events and trasitions
figure(1)
hold on
plot(diff(all_evts_r), 'k')
hline(peak_threshold, '--r')
plot(Rec_idx, 100, '*k')

% break them into recording sessions
t_start = Rec_idx(1:2:end-1);
t_end = Rec_idx(2:2:end);


plot([t_start ; t_end]', [50 50], '-b')

% convert identified peaks back into a time domain
Rec_intervals = all_evts_r(Rec_idx);
Inter_start = Rec_intervals(1:2:end-1);
Inter_end = Rec_intervals(2:2:end);

Rec_time_from_start = Rec_intervals  - csc.tvec(1);

for ii = 1: length(Inter_start)
    if ii  ==1
        time_since = minutes((Inter_start(ii)  - csc.tvec(1))/60);
    else
        time_since = minutes(((Inter_start(ii)- csc.tvec(1)) - (Inter_end(ii-1) - csc.tvec(1)))/60);
    end
    fprintf('\nDetected intervals %.0f started at %s, %s since last Ca2 recording, and was %.2f minutes long',ii, char(NLX_start + minutes((Inter_start(ii)  - csc.tvec(1))/60)),char(time_since),  (Inter_end(ii) - Inter_start(ii))/60)
end
fprintf('\n')

    fprintf('Number of recording sessions from Ca2+ MS file: %.0f\n', length(ms.timestamps))

    
%% Try to get times from imaging data
    
%     [~, peak_idx] = findpeaks(abs(diff(ms.time)), 

%% plot
figure(8)

% plot(csc.tvec(1:10000), csc.data(1,1:10000), csc.tvec(1:10000), csc.data(2,1:10000),evt.t{3}, '*k' )
% t_start = nearest_idx3(csc.tvec, evt.t{3}(1));
% t_end = nearest_idx3(csc.tvec, evt.t{3}(end));

%
plot(csc_r.tvec, csc_r.data(1,:))
% convert x axis to reasonable time units.  In this case hours
% x_val = get(gca, 'xtick');
% set(gca, 'xticklabel', round(((x_val - x_val(1))/60)/60,1))


hold on
plot(all_evts_r,max(csc_r.data(1,:)), '*k' )
% plot(evt.t{4},max(csc_r.data(1,:)), '*c' )


