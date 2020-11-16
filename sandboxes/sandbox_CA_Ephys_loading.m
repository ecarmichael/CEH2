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
    %     rmpath([PARAMS.code_base_dir '/io/neuralynx'])
    addpath(genpath('/Users/jericcarmichael/Documents/NLX_loaders_UNIX_2015'))
    disp('Version is greater than 2014b on UNIX so use updated loaders found here:')
    which Nlx2MatCSC
end

clear d os
%%  Laod some stuff
tic
load('ms.mat');
% load('SFP.mat');
% TS = MS_Load_TS('timestamp.dat');

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
cfg_ts = [];
% cfg_ts.fname = {'timestamp12.dat'};
cfg_ts.correct_time = 1;

TS= MS_Load_TS(cfg_ts);

all_TS_vals = [];
fprintf('\n****Comparing TS files to processed miniscope (ms) data\n')
for iT = 1:length(TS)
    if length(TS{iT}.system_clock{1}) == ms.timestamps(iT)
        disp([TS{iT}.filename   ':  ' num2str(length(TS{iT}.system_clock{1}))   ' - ms TS: ' num2str(ms.timestamps(iT))  '   ~ ' num2str(length(TS{iT}.system_clock{1}) / TS{iT}.cfg.Fs{1}) 's'])
    else
        warning(['TS do not match ms data' TS{iT}.filename   ':  ' num2str(length(TS{iT}.system_clock{1}))   ' - ms TS: ' num2str(ms.timestamps(iT))])
    end
    all_TS_vals = [all_TS_vals; TS{iT}.system_clock{end}];
end

% TS2 = MS_Load_TS('timestamp2.dat');

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

check = 1; check_fig =101;
tic


% load NLX events
cfg = [];
evt = LoadEvents(cfg);
evt.t{length(evt.t)+1} = unique(sort([evt.t{3} evt.t{4}]));
evt.label{length(evt.label)+1} = 'all_evt';

% load the NLX 'continuous sampled channels' "CSC"
cfg = [];
cfg.fc = {'CSC8.ncs'}%,'CSC6.ncs', 'CSC8.ncs'};
cfg.decimateByFactor = 8;
csc = LoadCSC(cfg);

% csc = restrict(csc, 1, 50000);

if check == 1
    len = 1:4000;
    figure(check_fig)
    ax(1) =  subplot(3,1,1);
    tvec = csc.tvec(len) - csc.tvec(1);
    offset = 0.0002;
    hold on
    for iC = 1:length(csc.label)
%         plot(tvec, csc.data(iC,len)+(offset*iC))cfg_filt = [];
cfg_filt.type = 'butter'; %Cheby1 is sharper than butter
cfg_filt.f  = [140 250]; % broad, could use 150-200?
cfg_filt.order = 4; %type filter order (fine for this f range)
cfg_filt.display_filter = 0; % use this to see the fvtool
csc_ripple = FilterLFP(cfg_filt, this_csc);
    end
end


% get the recording start time from NLX header
if isfield(csc.cfg.hdr{1}, 'TimeCreated')
    NLX_start = csc.cfg.hdr{1}.TimeCreated; % find the creation time as a string
    if contains(NLX_start, ':')
        NLX_start = duration(str2double(strsplit(NLX_start(end-8:end),':'))); % pull out hours:mins:sec and convert to a time
    else
        NLX_start = duration([NLX_start(end-5:end-4) ':' NLX_start(end-3:end-2) ':' NLX_start(end-1:end)]); % pull out hours:mins:sec and convert to a time
    end
end

toc


%% filters

% filter into the theta band
cfg_filt = [];
% cfg_filt.f = [5 11]; %setting theta (hertz)
cfg_filt.type = 'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt.f  = [1 5];
cfg_filt.order = 8; %type filter order
cfg_filt.display_filter = 1; % use this to see the fvtool
delta_csc = FilterLFP(cfg_filt, csc);

% filter into the theta band
cfg_filt = [];
% cfg_filt.f = [5 11]; %setting theta (hertz)
cfg_filt.type = 'cheby1'; %the type of filter I want to use via filterlfp
cfg_filt.f  = [6 11];
cfg_filt.order = 3; %type filter order
cfg_filt.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% reason.  Looks like trash with ord ~= 3 in cheby1. Butter and fdesign are also trash.
theta_csc = FilterLFP(cfg_filt, csc);

% add in the theta -  delta ratio (using just the filtered signals w/ 1s smoothing.
% t_d_ratio = smooth(abs(hilbert(theta_csc.data))./abs(hilbert(delta_csc.data)), theta_csc.cfg.hdr{1}.SamplingFrequency); % get the theta / delta ratio with some smoothing

% % alternative filtering (not used, just a check)
%
% Fs = csc.cfg.hdr{1}.SamplingFrequency;
% Wp = [ 6 11] * 2 / Fs;
% Ws = [ 4 13] * 2 / Fs;
% [N,Wn] = cheb1ord( Wp, Ws, 3, 20); % determine filter parameters
% [b_c1,a_c1] = cheby1(N,0.5,Wn); % builds filter
% csc_filtered = filtfilt(b_c1,a_c1,csc.data(1,:));



if check ==1
    figure(check_fig)
    ax(2) = subplot(3,1,2);
    
    hold on
    tvec = theta_csc.tvec(len) - theta_csc.tvec(1);
    hold on
    offset = 0.0002;
    for iC = 1:length(theta_csc.label)
%         plot(tvec, abs(hilbert(theta_csc.data(iC,len)))+(offset*iC))
        plot(tvec, (theta_csc.data(iC,len))+(offset*iC))
    end
    linkaxes(ax, 'x')
end



% identify major jumps in evts

%  all_jumps = diff(evt.t{5}) > (mean(diff(evt.t{5}) +0.5*std(diff(evt.t{5}))));
%  all_jumps(1) = 0; % correct for first jump;
%  jump_idx = find(all_jumps ==1);
%  rec_evt = [];
%  if sum(all_jumps) > 0 && sum(all_jumps) <2
%      fprintf('Jump found at time: %.0f\n', evt.t{5}(jump_idx))
%
%      rec_evt{1} = restrict(evt, evt.t{5}(1), evt.t{5}(jump_idx)); % add one index to compensate for the diff.
%      rec1_csc = restrict(csc, evt.t{5}(1), evt.t{5}(jump_idx));
%
%      rec_evt{2} = restrict(evt, evt.t{5}(jump_idx+1), evt.t{5}(end));
%      rec2_csc = restrict(csc, evt.t{5}(jump_idx+1), evt.t{5}(end));
%
%  elseif sum(all_jumps) >2
%
%      for iJ = length(jump_idx):-1:1
%          if iJ ==1
%              rec_evt{iJ} = restrict(evt, evt.t{5}(1), evt.t{5}(jump_idx(iJ)));
% %              rec_csc{iJ} = restrict(csc, evt.t{5}(1), evt.t{5}(jump_idx(iJ)));
%          else
%              rec_evt{iJ} = restrict(evt, evt.t{5}(jump_idx(iJ-1)), evt.t{5}(jump_idx(iJ)));
% %              rec_csc{iJ} = restrict(csc, evt.t{5}(jump_idx(iJ-1)), evt.t{5}(jump_idx(iJ)));
%          end
%      end
%  end

%% if the TSs align with the evt then add it in as a subfield [works for EVA only]
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

% maybe find jumps and fill them in?
TS_cam_mode = mode(diff(TS{iE}.system_clock{1}(2:end)));
% for iEvt = length(TS{iE}.system_clock{1}):-1:2

%     if TS{iE}.system_clock{1}(iEvt) - TS{iE}.system_clock{1}(iEvt-1)

% make a interpolated signal to see where things are missing
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for JISU DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % try plotting all of TS events as an array but with different colors for
% % each section
% all_ts = []; all_ys = []; this_val = 0.5;
% % colors = repmat({'r', 'b', 'g', 'c'},1,length(TS))
% for iT = 1:length(TS)
%     all_ts = [all_ts, TS{iT}.system_clock{1}'];
%     if this_val >0
%         all_ys = [all_ys, repmat(this_val,1, length(TS{iT}.system_clock{1}'))];
%         this_val = -.5;
%     elseif this_val <0
%                 all_ys = [all_ys, repmat(this_val,1, length(TS{iT}.system_clock{1}'))];
%                 this_val =.5 ;
%     end
% %     all_ys = [all_ys, repmat(iT, 1, length( TS{iT}.system_clock{1}'))];
% end

%% flag the longest part of the TS and flag it as a track sesson as it's own TS structure

for iT = 1:length(TS)
    all_TS_size(iT) = length(TS{iT}.system_clock{end});
    cam_ids(iT,:) = unique(TS{iT}.camera_id{end});
end

[largest_TS, idx] = max(all_TS_size);

fprintf('\nLargest TS segment is %.0f samples ~ %0.2f mins  cam ID: %.0f\n', largest_TS, (largest_TS/TS{idx}.cfg.Fs{end})/60, cam_ids(idx))

% remove the largest from the TS structure and put it in a new one

% make new TS for track
TS_track{1} = TS{idx};
% now remove it
TS(idx) = [];


% try to segment the ms structure
ms_seg = MS_segment_ms_sandbox(ms);


% remove the track segment
ms_seg = MS_remove_data_sandbox(ms_seg, [idx]);
%% identify peaks in  diff(evt.t{5}) marking transitions in the camera TTLs
peak_threshold =  (mean(diff(evt.t{5}) +0.05*std(diff(evt.t{5}))));
min_dist = 10;
[Rec_peak, Rec_idx] = findpeaks(diff(evt.t{3}), 'minpeakheight',peak_threshold, 'minpeakdistance', min_dist);
fprintf(['\nDetected %.0f trigger transitions treating this as %.0f distinct recordings\n'], length(Rec_idx), length(Rec_idx))


% plot the diff and the detected peaks as a check.
if check == 1
    figure(check_fig)
    hold on
    subplot(3,1,3)
    plot(diff(evt.t{3}), 'k')
    hline(peak_threshold, '--r')
    % plot(Rec_idx, 100, '*k')
    for iRec = 1:length(Rec_idx)
        % text(Rec_idx(iRec),Rec_peak(iRec),num2str(iRec))
        text(Rec_idx(iRec),peak_threshold,num2str(iRec))
        
    end
end


% t_start = Rec_idx(1:2:end-1);
% t_end = Rec_idx(2:2:end);
% plot([t_start ; t_end]', [50 50], '-b')


%% check the transitions for jumps and compare them to the TS times.  It seems like the TTL can have gitter around the transition periods.
t_idx = 3;

for iRec = 1:length(Rec_idx)
    
    low_val = 10;
    high_val = -10;
    idx_low = NaN; % initialize with something.
    idx_high = NaN;
    
    if iRec == 1
        low_val = 0;  % only move the low_val 'forward' from zero
        while ~isempty(idx_low)
            temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-low_val), evt.t{t_idx}(Rec_idx(iRec+1)));
            idx_low =find(diff(temp_evt.t{t_idx}(1:20)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            %     idx_high =find(diff(temp_evt.t{t_idx}(end-50:end)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            low_val = low_val - 1;
        end
        
        
        while ~isempty(idx_high)
            temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)), evt.t{t_idx}(Rec_idx(iRec+1)-high_val));
            %     idx_low =find(diff(temp_evt.t{t_idx}(1:50)) > mode(diff(temp_evt.t{t_idx}))*1.05)
            idx_high =find(diff(temp_evt.t{t_idx}(end-20:end)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            high_val = high_val +1;
        end
        
        disp(['Corrected indexting ' num2str(0) ' - ' num2str(high_val)])
        
        
        
    elseif iRec == length(Rec_idx)
        high_val = 0; % only move backwards from the end.
        
        while ~isempty(idx_low)
            temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-low_val), evt.t{t_idx}(end));
            idx_low =find(diff(temp_evt.t{t_idx}(1:20)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            %     idx_high =find(diff(temp_evt.t{t_idx}(end-50:end)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            low_val = low_val - 1;
        end
        
        while ~isempty(idx_high)
            temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)), evt.t{t_idx}(end)-high_val);
            %     idx_low =find(diff(temp_evt.t{t_idx}(1:50)) > mode(diff(temp_evt.t{t_idx}))*1.05)
            idx_high =find(diff(temp_evt.t{t_idx}(end-20:end)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            high_val = high_val +1;
        end
        
        
        disp(['Corrected indexting ' num2str(low_val) ' - ' num2str(0)])
        
    else
        
        % loop until you find the index offset that gives no jump in time. for
        % start
        while ~isempty(idx_low)
            temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-low_val), evt.t{t_idx}(Rec_idx(iRec+1)));
            idx_low =find(diff(temp_evt.t{t_idx}(1:15)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            %     idx_high =find(diff(temp_evt.t{t_idx}(end-50:end)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            low_val = low_val - 1;
        end
        
        % loop until you find the index offset that gives no jump in time. for
        % end
        while ~isempty(idx_high)
            temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-low_val), evt.t{t_idx}(Rec_idx(iRec+1)-high_val));
            %     idx_low =find(diff(temp_evt.t{t_idx}(1:50)) > mode(diff(temp_evt.t{t_idx}))*1.05)
            idx_high =find(diff(temp_evt.t{t_idx}(end-15:end)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            high_val = high_val +1;
        end
        
        
        disp(['Corrected indexting ' num2str(low_val) ' - ' num2str(high_val)])
        
    end
    
    times_to_use(iRec,:) = [low_val+1, high_val-1]; % keep the good values for restrictions in the next cell.
    
    if iRec == length(Rec_idx)
        temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-low_val+1), evt.t{t_idx}(end)-high_val-1);
    else
        temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-low_val+1), evt.t{t_idx}(Rec_idx(iRec+1)-high_val-1));
    end
    biggest_jumps(iRec) = max(diff(temp_evt.t{t_idx}));
end

for iRec = 1:length(Rec_idx)
    fprintf('Evt: %.0f Best offsets: start = %.0f  end = %.0f largest jump = %.3f sec\n', iRec, times_to_use(iRec,1), times_to_use(iRec, 2), biggest_jumps(iRec))
    
end


%% make some EVT blocks corresponding to the transitions and chop the data

% allocate nCells
rec_evt = cell(length(Rec_idx),1);
rec_csc = cell(length(Rec_idx),1);
rec_theta = cell(length(Rec_idx),1);

% restrict to the recording period and put them in cells. for events:
% 'evt', raw lfp : 'csc', filtered lfp = 'theta'
for iRec = 1:length(Rec_idx)
    if iRec < length(Rec_idx)
        rec_evt{iRec} = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-times_to_use(iRec,1)), evt.t{t_idx}(Rec_idx(iRec+1)-times_to_use(iRec,2))); % restrict the NLX evt struct to ms TTL periods
        rec_csc{iRec} = restrict(csc, evt.t{t_idx}(Rec_idx(iRec)+2), evt.t{t_idx}(Rec_idx(iRec+1)-1)); % same for the csc
        % same for filtered traces (if they exist)
        if exist('delta_csc', 'var')
            rec_delta{iRec} = restrict(delta_csc, evt.t{t_idx}(Rec_idx(iRec)+2), evt.t{t_idx}(Rec_idx(iRec+1)-1)); % same for the csc
        end
        if exist('theta_csc', 'var')
            rec_theta{iRec} = restrict(theta_csc, evt.t{t_idx}(Rec_idx(iRec)+2), evt.t{t_idx}(Rec_idx(iRec+1)-1)); % same for the csc
        end
    else
        rec_evt{iRec} = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-times_to_use(iRec,1)), evt.t{t_idx}(end)); % restrict the NLX evt file (last only)
        rec_csc{iRec} = restrict(csc, evt.t{t_idx}(Rec_idx(iRec)+2), evt.t{t_idx}(end)); % same for csc
        
        % for filtered traces.
        if exist('delta_csc', 'var')
            rec_delta{iRec} = restrict(delta_csc, evt.t{t_idx}(Rec_idx(iRec)+2), evt.t{t_idx}(end)); % same for csc
        end
        if exist('theta_csc', 'var')
            rec_theta{iRec} = restrict(theta_csc, evt.t{t_idx}(Rec_idx(iRec)+2), evt.t{t_idx}(end)); % same for csc
        end
        
        
    end
end

%% check for jumps in TS files
for iT = 1:length(TS)
    fprintf('TS %s mode diff = %.0f max diff = %.0f\n', TS{iT}.filename, mode(diff(TS{iT}.system_clock{end})), max(diff(TS{iT}.system_clock{end})))
end

%% check length of TSs
all_evt = 0;
for iRec = 1:length(rec_evt)
    %     disp(['Rec ' num2str(iRec)])
    for this_evt = length(rec_evt{iRec}.label) % correct for start and stop recording.
        fprintf('Number of evts id: %.0f  =   %.0f samples at %0.2f Hz. Start at: %s for ~%.1fs \n',iRec, length(rec_evt{iRec}.t{this_evt}),1/(median(diff(rec_evt{iRec}.t{this_evt}))),...
            char(NLX_start + minutes((rec_evt{iRec}.t{this_evt}(1)  - csc.tvec(1))/60)), length(rec_evt{iRec}.t{this_evt})/(1/(median(diff(rec_evt{iRec}.t{this_evt})))))
    end
    all_evt = all_evt + length(rec_evt{iRec}.t{this_evt});
    all_evt_lens(iRec) = length(rec_evt{iRec}.t{this_evt});
end

all_TS = 0;
for iRec = 1:length(TS)
    %     disp(['TS ' num2str(iRec)])
    fprintf('Number of Scope TS id: %.0f  =   %.0f  at %0.2fHz for %.1f sec\n',iRec, length(TS{iRec}.system_clock{1}), 1/(median(diff(TS{iRec}.system_clock{1}(2:end)))*0.001),...
        length(TS{iRec}.system_clock{1})/ (1/(median(diff(TS{iRec}.system_clock{1}(2:end)))*0.001)))
    all_TS = all_TS + length(TS{iRec}.system_clock{1});
    all_TS_len(iRec) = length(TS{iRec}.system_clock{1});
    
end

fprintf('All EVT: %.0f  All TS: %.0f', all_evt, all_TS)

disp('All evt')
for this_evt = 3:length(evt.label) % correct for start and stop recording.
    fprintf('Number of evt evts id: %.0f  =   %.0f at %0.2f Hz\n',this_evt, length(evt.t{this_evt}),1/(median(diff(evt.t{this_evt}))))
end

% fprintf('Number of NLX events: %.0f, Number of Scope TS: %.0f, Difference: %.0f\n', length(evt.t{this_evt}), length(TS.SYSCLOCK{this_cam}(2:end)), length(evt.t{this_evt}) -  length(TS.SYSCLOCK{this_cam}(2:end)));

% compare
disp('Compare')

for iRec = 1:length(rec_evt)
    %     disp(['Rec ' num2str(iRec)])
    for this_evt = length(rec_evt{iRec}.label) % correct for start and stop recording.
        fprintf('Evts id: %.0f = %.0f samples fs ~ %.1f time: %0.2f sc || TS id: %.0f = %.0f samples fs ~ %.1f time: %0.2f sc\n',iRec, length(rec_evt{iRec}.t{this_evt}),mode(diff(rec_evt{iRec}.t{this_evt}))*1000,length(rec_evt{iRec}.t{this_evt})/(1/(median(diff(rec_evt{iRec}.t{this_evt})))),...
            iRec, length(TS{iRec}.system_clock{end}), mode(diff(TS{iRec}.system_clock{end})),...
            length(TS{iRec}.system_clock{1})/ (1/(median(diff(TS{iRec}.system_clock{1}(2:end)))*0.001)))
    end
    evt_TS_diff(iRec) = length(rec_evt{iRec}.t{this_evt}) - length(TS{iRec}.system_clock{end});
end

evt_TS_diff % print the offset
%% Add in csc data to the ms data strucuture.  


% append restricted csc files
ms_seg = MS_append_data_sandbox(ms_seg, 'csc', rec_csc');

% appened a theta filtered signal

ms_seg = MS_append_data_sandbox(ms_seg, 'theta_csc', rec_theta');

ms_seg = MS_append_data_sandbox(ms_seg, 'delta_csc', rec_delta');

%% Plot some examples of segments

Chans = 1:5; % channels to plot
c_ord = linspecer(length(Chans)); % nice colours.
plot_type = '3d'; % ploting style can be '2d' or '3d'

for iRec = 1:3
    figure(iRec)
    ax(1) =subplot(2,1,1);
    timein = (ms_seg.theta_csc{iRec}.tvec - ms_seg.theta_csc{iRec}.tvec(1)); % just to fix the timing offset between them back to ebing relative to this segment.
    % timein = timein
    
    plot(timein, abs(hilbert(ms_seg.theta_csc{iRec}.data)), '--r');
    hold on
    plot(timein,ms_seg.theta_csc{iRec}.data, '-b' );
    xlim([timein(1), timein(end)])
    
    ax(2) =subplot(2,1,2);
    hold on
    time_in2 = ms_seg.time{iRec} - ms_seg.time{iRec}(1);
    switch plot_type
        % 2d
        case '2d'
            for iC = 1:length(Chans)
                plot(time_in2*0.001, ms_seg.RawTraces{iRec}(:,iC), 'color', c_ord(iC,:))
            end
            % 3d
        case '3d'
            for iC = 1:length(Chans)
                plot3(time_in2*0.001, repmat(iC,size(ms_seg.RawTraces{iRec},1),1), ms_seg.RawTraces{iRec}(:,iC), 'color', c_ord(iC,:))
            end
            view([0 45])
    end
    
    xlim([time_in2(1)*0.001 time_in2(end)*0.001])
    linkaxes(ax, 'x')
    ax = [];
end


%% try getting a spectrogram for sleep state detection (if this works move it to before the data is segmented)
cfg_spec.win = 512*2;
cfg_spec.noverlap = floor(cfg_spec.win/10);



[~,F,T,P] = spectrogram(csc.data(1:round(length(csc.data)/2)),hamming(cfg_spec.win),cfg_spec.noverlap,1:0.1:40,csc.cfg.hdr{1}.SamplingFrequency, 'power');

% chronux
movingwin=[1 0.5]; % set the moving
params = [];
params.Fs = csc.cfg.hdr{1}.SamplingFrequency; % sampling frequency
params.fpass=[1 100]; % frequencies of
params.tapers=[3 5]; % tapers
% params.trialave=1; % average over trials
% params.err=0; % no error
[S,t,f]=mtspecgramc(csc.data(1:round(length(csc.data)/2)),movingwin,params);

figure(123)
   imagesc(t,f,10*log10(S)');axis xy; colorbar; title('Spectrogram');
% axis xy
% caxis([-130 -70])


figure(122)
ax1 = imagesc(T/60/60,F,10*log10(P)); % converting to dB as usual
set(ax1, 'AlphaData', ~isinf(10*log10(P)))
%         set(gca,'FontSize',28);
axis xy; xlabel('Time (hr)'); ylabel('Frequency (Hz)');
% xlim([10 12])
caxis([-130 -70])
hold on
plot(csc.tvec(1:round(length(csc.data)/2)) - csc.tvec(1), (csc.data(1:round(length(csc.data)/2))*10000)+(max(F)/2), 'k')
plot(csc.tvec(1:round(length(csc.data)/2)) - csc.tvec(1), (theta_csc.data(1:round(length(csc.data)/2))*20000)+10, 'c')
plot(csc.tvec(1:round(length(csc.data)/2)) - csc.tvec(1), (delta_csc.data(1:round(length(csc.data)/2))*20000)+3, 'm')
% plot(csc.tvec(1:round(length(csc.data)/2)) - csc.tvec(1), (t_d_ratio(1:round(length(csc.data)/2))*0.5), '--r')


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


