function [data_out] = Master_OE_Preprocess(data_dir, type)
%% Master_OE_Preprocess: Loads and preprocesses the Radial/NOL maze data recorded using a combination of OE ephys with miniscope video tracking
%
%
%
%    Inputs: 
%    - data_dir: [path]   location of the dir containing all of the OE
%    recordings and the miniscope 'day' folder containing the encoding and
%    recall. 
%    
%    - type: [string]  
%
%    Outputs: 
%    -
%
%
%
%
% EC 2024-10-21   initial version 
%
%
%
%% initialize
p = strsplit(data_dir, filesep); 
if contains(p{end}, 'Rad') || contains(p{end}, 'AA')
    type = 'Radial';
elseif contains(p{end}, 'NOL')
    type = 'NOL';
end
    
fprintf('%s: Type is set to %s...\n', mfilename, type)

cd(data_dir)

f_list = dir(fullfile(data_dir, '**\spike_times.npy')); 
kilo_dir = f_list.folder; 

f_list = dir(fullfile(data_dir, '**\*A.events')); 
OE_dir = f_list.folder; 

% where are the tracking files

f_list = dir(fullfile(data_dir, '**\timestamps.csv'));
DLC_dir = []; 

for ii = 1:length(f_list)
    if contains(f_list(ii).folder, 'minicam1')
        DLC_dir{end+1} = f_list(ii).folder; 
    end
end


if isunix
    
else
    save_dir = ['C:\Users\' getenv('username') '\Williams Lab Dropbox\Williams Lab Team Folder\Eric\' type '\inter\'];
end
%% Load and preprocess the data



%% get the meta or create it. 
cd(OE_dir)
[p_dir] = fileparts(cd); 
cd(p_dir)
MS_Write_meta_RAD
meta = MS_Load_meta;

%% load the behaviour information. 
cd(p_dir)
rad_name = dir('Radial*.m');

if ~isempty(rad_name) && length(rad_name) <2
    run(rad_name.name)
        
    target = Rad.(['D' rad_name.name(end-11:end-2)]).correct;
    trl = Rad.(['D' rad_name.name(end-11:end-2)]).(rad_name.name(8:end-13));
    
    % Encoding trials
    
    Enc_iv = iv(trl.encode.tstart(1:4), trl.encode.tend(1:4));
    
    iti_s = [trl.encode.tstart(1:4)-60, trl.encode.tstart(5)];
    iti_e = [trl.encode.tstart(1:4), trl.encode.tend(5)];
    
    Enc_iti_iv = iv(iti_s, iti_e);
    
    % Recall trials
    
    Rec_iv = iv(trl.recall.tstart(1:4), trl.recall.tend(1:4));
    
    iti_s = [trl.recall.tstart(1:4)-60, trl.recall.tstart(5)];
    iti_e = [trl.recall.tstart(1:4), trl.recall.tend(5)];
    
    Rec_iti_iv = iv(iti_s, iti_e);
elseif isempty(rad_name) 
    error('No Radial* behaviour log')
elseif length(rad_name) > 1
    error('More than one Radial behaviour log')
    
end
    

%% load the event
cd(OE_dir)
mess = OE_LoadMessages; 

evts = OE_LoadEvents; 

if sum(contains(evts.label, meta.TTL_rec)) == 1
    evts.label{(contains(evts.label,  meta.TTL_miniscope_frame))} = 'TTL_frame'; 
    evts.label{(contains(evts.label,  meta.TTL_rec))} = 'TTL_rec'; 
    
    % add lasers if present. 
    if sum(contains(evts.label,  meta.TTL_LED)) > 0
        evts.label{(contains(evts.label,  meta.TTL_LED))} = 'TTL_LED';
        LED_idx = find(contains(evts.label, 'TTL_LED'));
        
        % convert laser TTLs into blocks
        keep_idx = diff(evts.t{LED_idx}) > 4;
        keep_idx = [1 ; find(keep_idx(1:end))+1]';
        
        
        laser_on = evts.t{LED_idx}(keep_idx);
        laser_off = evts.t{LED_idx}([keep_idx(2:end)-1 length(evts.t{LED_idx})]);
        
        LED_IV = iv(laser_on, laser_off);
    end
    
    
    frame_idx = find(contains(evts.label, 'TTL_frame'));
    trig_idx = find(contains(evts.label, 'TTL_rec'));
    


    for ii = 1:length(evts.label)
        dur(ii) = evts.t{ii}(end) - evts.t{ii}(1);
    end
    
    
    [~, idx] = findpeaks(diff(evts.t{frame_idx}), 'MinPeakHeight', 1000);
    
    enc_t = [evts.t{frame_idx}(1) evts.t{frame_idx}(idx)]; 
    rec_t = [evts.t{frame_idx}(idx+1) evts.t{frame_idx}(end)]; 

    [~, max_TTL] = max(dur); 
    start_t = evts.t{max_TTL}(1);
    end_t = evts.t{max_TTL}(end); 
    
    % get the sleep times 
    sleep_idx = find(contains(mess.label, 'sleep')); 
    sleep_end_idx = find(contains(mess.label, 'sleep end')); 
    
    if length(sleep_idx) >1
        sleep_idx = sleep_idx(1); 
    end
    
    if sleep_idx == sleep_end_idx % if only sleep end was logged;
        sleep_end_t = mess.t{sleep_end_idx}; 
        sleep_t = sleep_end_t - 60*4*60; % subtract 4 hours in seconds. 
        sleep_IV = iv(sleep_t, sleep_end_t); 
    elseif ~isempty(sleep_idx) && isempty(sleep_end_idx) % case where sleep was logged but end was not. 
        sleep_t = mess.t{sleep_idx}; 
        sleep_end_t = sleep_t + 60*4*60; % subtract 4 hours in seconds. 
        sleep_IV = iv(sleep_t, sleep_end_t); 
    else % case where both sleep and sleep end were logged. 
        sleep_IV = iv(mess.t{sleep_idx}, mess.t{sleep_end_idx}); 
    end


    
else
    error('No trigger TTL found. Required for synchronizing data')
end


fprintf('<strong>%s</strong>: OE frame TTL (<strong>''%s''</strong>), duration = %.0fsec (%2.1fhr)\n', mfilename,evts.label{max_TTL}, end_t - start_t, (end_t - start_t)/60/60)
fprintf('<strong>%s</strong>: Encode duration = %.0fsec (%2.1fmin)\n', mfilename,diff(enc_t), diff(enc_t)/60);
fprintf('<strong>%s</strong>: Sleep duration = %.0fsec (%2.1fmin = %2.1fhrs)\n', mfilename,sleep_IV.tend - sleep_IV.tstart, (sleep_IV.tend - sleep_IV.tstart)/60, ((sleep_IV.tend - sleep_IV.tstart)/60)/60)
fprintf('<strong>%s</strong>: Recall duration = %.0fsec (%2.1fmin)\n', mfilename,diff(rec_t), diff(rec_t)/60)


%% load spikes if they are there
cd(kilo_dir)

S = OE_phy2TS;

% convert to time. 
for ii = 1:length(S.t)
   S.t{ii} = S.t{ii}./30000; 
end


%% load the CSC data. 
cd(OE_dir)

cfg_csc = [];
cfg_csc.fc =  {meta.goodCSC};%cfg.csc_chan;
cfg_csc.desired_sampling_frequency = 2000;
csc = OE_old_csc2TSD(cfg_csc);

csc.cfg.hdr{1} = rmfield(csc.cfg.hdr{1}, {'ts', 'nsamples', 'recNum'});

if ~strcmp(meta.EMG, 'NA')
    cfg_csc = [];
    cfg_csc.fc =  {meta.OE_EMG};
    cfg_csc.desired_sampling_frequency = 2000;
    emg = OE_old_csc2TSD(cfg_csc);
    
    
    %
    csc.data(end+1,:) = emg.data;
    csc.label{end+1} = 'EMG';
    csc.cfg.hdr{end+1} = emg.cfg.hdr{1};
end

%% check for accel

if isfield(meta, 'OE_acc')

cfg_acc = [];
cfg_acc.fc =  meta.OE_acc;%cfg.csc_chan;
cfg_acc.desired_sampling_frequency = 2000;
acc = OE_old_csc2TSD(cfg_acc);
% fs = cfg_acc.desired_sampling_frequency; 
% eulZYX_df_t(1,:) = conv2(diff(acc.data(1,:)),gausswin(fs, 3),'same'); 
% eulZYX_df_t(2,:) = conv2(diff(acc.data(2,:)),gausswin(fs, 3),'same');
% eulZYX_df_t(3,:) = conv2(diff(acc.data(3,:)),gausswin(fs, 3),'same');
% 
% eulZYX_df = [eulZYX_df_t, eulZYX_df_t(:,end)];

acc.data(4,:) = movmean(sqrt((acc.data(1,:).^2) + (acc.data(2,:).^2) + (acc.data(3,:).^2)),cfg_acc.desired_sampling_frequency/4);  %   sqrt(movmean(sum(abs(acc.data)).^ 2, cfg_acc.desired_sampling_frequency/4)); 
acc.label = {'AUX1', 'AUX2', 'AUX3', 'RMS'};
acc.cfg.hdr{1} = rmfield(acc.cfg.hdr{1}, {'ts', 'nsamples', 'recNum'});

end

%%  Get the position data if it is there

if ~isempty(DLC_dir)
    
    if length(DLC_dir) <2
    
        [pos, ~] = MS_DLC2TSD(DLC_dir, [], meta.conv_fact);
        if pos.tvec(1) <0
            pos.tvec = pos.tvec + abs(pos.tvec(1));
        elseif pos.tvec(1) >0
            pos.tvec  = pos.tvec - abs(pos.tvec(1));
        end
    
    else
        this_pos = []; 
        all_tvec = []; all_data = [];
        for ii = 1:length(DLC_dir)
            [this_pos, ~] = MS_DLC2TSD(DLC_dir{ii}, [],meta.conv_fact);
            if this_pos.tvec(1) <0
                this_pos.tvec = this_pos.tvec + abs(this_pos.tvec(1));
            elseif this_pos{ii}.tvec(1) >0
                this_pos.tvec  = this_pos.tvec - abs(this_pos.tvec(1));
            end
            
            if ii == 1
                enc.pos = this_pos;
                enc.pos.tvec = linspace(enc_t(1), enc_t(2), length(this_pos.tvec));
                enc.pos.label = {'x', 'y', enc.pos.label{2:end}};
            elseif ii == 2
                rec.pos = this_pos;
                rec.pos.tvec = linspace(rec_t(1), rec_t(2), length(this_pos.tvec));
                rec.pos.label = {'x', 'y', rec.pos.label{2:end}};

            end
%             all_tvec = [all_tvec, this_pos{ii}.tvec']; 
%             all_data = [all_data, this_pos{ii}.data];
%             this_pos{1} = tsd(all_tvec', all_data, 'label',this_pos{1}.label);
%          this_pos{1}.units = this_pos{1}.units; 

        end
        
        % merge the pos files
%         pos = tsd(all_tvec, all_data, 'label',this_pos{1}.label);
%         pos.units = this_pos{1}.units; 
        
    end
end
fprintf('<strong>%s</strong>: Encode Pos duration = %.0fs = %.1fmin\n', mfilename, enc.pos.tvec(end) - enc.pos.tvec(1),(enc.pos.tvec(end) - enc.pos.tvec(1))/60)
fprintf('<strong>%s</strong>: Recall Pos duration = %.0fs = %.1fmin\n', mfilename, rec.pos.tvec(end) - rec.pos.tvec(1), (rec.pos.tvec(end) - rec.pos.tvec(1))/60)

pos = enc.pos; 
pos.tvec = [enc.pos.tvec, rec.pos.tvec]; 
pos.data = [enc.pos.data, rec.pos.data]; 
clear this_pos

% speed = getLinSpd([],pos); % linear speed

%% align the behav times to OE times
% 
% fprintf('Pos dur(%.2fsec) - TTLs (%0.2fsec), diff: %0.3f\n', pos.tvec(end) - pos.tvec(1), end_t - start_t, (pos.tvec(end) - pos.tvec(1)) -  (end_t - start_t)) ; 
% 
% if ((pos.tvec(end) - pos.tvec(1)) -  (end_t - start_t)) >1
%     error('position and TTLs differ by more than 1sec')
% end
% 
% pos.tvec = pos.tvec + start_t; 

%% restrict to encoding, sleep, recall
enc.csc = restrict(csc, enc_t(1), enc_t(2)); 
sleep.csc = restrict(csc, enc_t(2), rec_t(1)); 
rec.csc = restrict(csc, rec_t(1), rec_t(2)); 

if exist('emg', 'var')

enc.emg = restrict(emg, enc_t(1), enc_t(2)); 
sleep.emg = restrict(emg, enc_t(2), rec_t(1)); 
rec.emg = restrict(emg, rec_t(1), rec_t(2)); 
end

if exist('acc', 'var')
enc.acc = restrict(acc, enc_t(1), enc_t(2)); 
sleep.acc = restrict(acc, enc_t(2), rec_t(1)); 
rec.acc = restrict(acc, rec_t(1), rec_t(2)); 
end
% 
enc.S = restrict(S, enc_t(1), enc_t(2)); 
sleep.S = restrict(S, enc_t(2), rec_t(1)); 
rec.S = restrict(S, rec_t(1), rec_t(2)); 

% update the trial information
enc.TRL_iv = Enc_iv; 
enc.TRL_iv.tstart = Enc_iv.tstart+ enc.csc.tvec(1); 
enc.TRL_iv.tend = Enc_iv.tend+ enc.csc.tvec(1); 

enc.ITI_iv = Enc_iti_iv; 
enc.ITI_iv.tstart = Enc_iti_iv.tstart+ enc.csc.tvec(1); 
enc.ITI_iv.tend = Enc_iti_iv.tend+ enc.csc.tvec(1); 

rec.TRL_iv = Rec_iv; 
rec.TRL_iv.tstart = Rec_iv.tstart+ rec.csc.tvec(1); 
rec.TRL_iv.tend = Rec_iv.tend+ rec.csc.tvec(1); 

rec.ITI_iv = Rec_iti_iv; 
rec.ITI_iv.tstart = Rec_iti_iv.tstart+ rec.csc.tvec(1); 
rec.ITI_iv.tend = Rec_iti_iv.tend+ rec.csc.tvec(1); 

%% check the tracking and trial intervals
figure(99)
clf
subplot(3,4,1:2)
hold on
plot(enc.pos.tvec, enc.pos.data(1:2,:), 'k')
plot(enc.acc.tvec, enc.acc.data(4,:)*500, 'r')

vline(enc.TRL_iv.tstart, '.-r');
vline(enc.TRL_iv.tend, '.-k');
vline(enc.ITI_iv.tstart, '--g'); 
vline(enc.ITI_iv.tend, '--m');
vline(enc_t, '--b');


subplot(3,4,3:4)
hold on
plot(rec.pos.tvec, rec.pos.data(1:2,:), 'k')
plot(rec.acc.tvec, rec.acc.data(4,:)*500, 'r')

vline(rec.TRL_iv.tstart, '.-r');
vline(rec.TRL_iv.tend, '.-k');
vline(rec.ITI_iv.tstart, '--g'); 
vline(rec.ITI_iv.tend, '--m');
vline(rec_t, '--b');

for ii = 1:4
    this_e = restrict(enc.pos, enc.TRL_iv.tstart(ii), enc.TRL_iv.tend(ii)); 
    this_e_iti = restrict(enc.pos, enc.ITI_iv.tstart(ii), enc.ITI_iv.tend(ii)); 
    
    this_r = restrict(rec.pos, rec.TRL_iv.tstart(ii), rec.TRL_iv.tend(ii));
    this_r_iti = restrict(rec.pos, rec.ITI_iv.tstart(ii), rec.ITI_iv.tend(ii)); 
    
    subplot(3,4,ii+4)
    hold on
    plot(this_e.data(1,:), this_e.data(2,:), '.r')
    plot(this_e_iti.data(1,:), this_e_iti.data(2,:), '.k')

    subplot(3,4,ii+8)
    hold on
    plot(this_r.data(1,:), this_r.data(2,:), '.r')
    plot(this_r_iti.data(1,:), this_r_iti.data(2,:), '.k')
    
    if ii ==1 
        legend({'TRL', 'ITI'})
    end
end
%% generate the hypnogram

csc_temp = sleep.csc;
csc_temp.data = csc_temp.data(1,:); 
csc_temp.label(2:end) = [];

csc_temp.cfg.hdr(2:end) = [];

if exist('acc', 'var')

    acc_temp = sleep.acc; 
    acc_temp.data = acc_temp.data(4,:); 
    acc_temp.label(2:end) = [];
    acc_temp.label{1} = sleep.acc.label{4}; 

hypno = dSub_Sleep_screener(1, csc_temp, acc_temp, []); 

else

hypno = dSub_Sleep_screener([], csc_temp, sleep.emg, []); 
end
close(221)

clear csc_temp acc_temp
%% save it all for output

out = [];
out.meta = meta; 
% out.S = S;
out.Encode = enc;
out.Sleep = sleep;
out.Recall = rec; 
out.evts = evts;
out.enc_t = enc_t;
out.rec_t = rec_t; 
% out.csc = csc;
% out.pos = pos;
out.hypno = hypno;
out.history.function{1} = mfilename;
out.history.date{1} = date;
out.history.cfg{1} = [];  


fname = [meta.subject '_' strrep(meta.date, '-', '_') '_' meta.session]; 

save([save_dir filesep fname], 'out')



