function out = Master_NOL_preprocess(cfg_in, Kilo_dir, OE_dir, DLC_dir, save_dir)
%% MS_preprocess_NOL:  Extracts the position, spike, and CSC data for an OE session.
%
%
%
%    Inputs:
%    - cfg_in: [struct]  configuration parameters to overwrite the defaults
%
%    - kilo_dir: [path]  directory with the processed kilosort/phy2
%    output
%
%    -  oE_dir: [path]  directory with the OE continuous files
%
%    - DLC_dir: [path] directory/directories with the DLC outputs. Can be a
%    cell array to process multiple directories for one session.
%
%    Outputs:
%    -
%
%
%
%
% EC 2023-06-29   initial version
%
%
%
%% initialize

cfg_def = [];
cfg_def.conv_fac = [340/38  350/38];
cfg_def.behav_label = 'Body';
cfg_def.TTL_trig = '6'; 
cfg_def.TTL_frame = '4'; 
cfg_def.csc_chan = {'CH44'};
cfg_def.emg_chan = {'CH18'}; 
cfg_def.acc_chan = {'CH65', 'CH66', 'CH67'};


cfg = ProcessConfig(cfg_def, cfg_in);


CSC = [];
S = [];
pos = [];
evt = [];
acc = [];

%% get the meta or create it. 
cd(OE_dir)
[p_dir] = fileparts(cd); 
cd(p_dir)
MS_Write_meta_NOL
meta = MS_Load_meta;


%% load the event
cd(OE_dir)
mess = OE_LoadMessages; 

evts = OE_LoadEvents; 

if sum(contains(evts.label, cfg.TTL_trig)) == 1
    evts.label{(contains(evts.label, cfg.TTL_frame))} = 'TTL_frame'; 
    evts.label{(contains(evts.label, cfg.TTL_trig))} = 'TTL_trig'; 
    
    frame_idx = find(contains(evts.label, 'TTL_frame'));
    trig_idx = find(contains(evts.label, 'TTL_trig'));

    for ii = 1:length(evts.label)
        dur(ii) = evts.t{ii}(end) - evts.t{ii}(1); 
    end
    
    [~, idx] = findpeaks(diff(evts.t{frame_idx}), 'MinPeakHeight', 1000);
    
    enc_t = [evts.t{frame_idx}(1) evts.t{frame_idx}(idx)]; 
    rec_t = [evts.t{frame_idx}(idx+1) evts.t{frame_idx}(end)]; 

    [~, max_TTL] = max(dur); 
    start_t = evts.t{max_TTL}(1);
    end_t = evts.t{max_TTL}(end); 
    
else
    error('No trigger TTL found. Required for synchronizing data')
end


fprintf('<strong>%s</strong>: OE frame TTL (<strong>''%s''</strong>), duration = %.0fsec (%2.1fhr)\n', mfilename,evts.label{max_TTL}, end_t - start_t, (end_t - start_t)/60/60)
fprintf('<strong>%s</strong>: Encode duration = %.0fsec (%2.1fmin)\n', mfilename,diff(enc_t), diff(enc_t)/60)
fprintf('<strong>%s</strong>: Recall duration = %.0fsec (%2.1fmin)\n', mfilename,diff(rec_t), diff(rec_t)/60)


%% load spikes if they are there
cd(Kilo_dir)

S = OE_phy2TS;

% convert to time. 
for ii = 1:length(S.t)
   S.t{ii} = S.t{ii}./30000; 
end


%% load the CSC data. 
cd(OE_dir)

cfg_csc = [];
cfg_csc.fc =  {meta.OE_goodCSC};%cfg.csc_chan;
cfg_csc.desired_sampling_frequency = 2000;
csc = OE_old_csc2TSD(cfg_csc);


cfg_csc = [];
cfg_csc.fc =  {meta.OE_EMG};
cfg_csc.desired_sampling_frequency = 2000;
emg = OE_old_csc2TSD(cfg_csc);

% 
csc.data(end+1,:) = emg.data;
csc.label{end+1} = 'EMG';
csc.cfg.hdr{end+1} = emg.cfg.hdr{1}; 

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

acc.data(4,:) = sqrt(movmean(sum(abs(acc.data)).^ 2, cfg_acc.desired_sampling_frequency/4)); 

end

%%  Get the position data if it is there

if ~isempty(DLC_dir)
    
    if length(DLC_dir) <2
    
        [pos, ~] = MS_DLC2TSD(DLC_dir, [], cfg.conv_fac);
        if pos.tvec(1) <0
            pos.tvec = pos.tvec + abs(pos.tvec(1));
        elseif pos.tvec(1) >0
            pos.tvec  = pos.tvec - abs(pos.tvec(1));
        end
    
    else
        this_pos = []; 
        all_tvec = []; all_data = [];
        for ii = 1:length(DLC_dir)
            [this_pos, ~] = MS_DLC2TSD(DLC_dir{ii}, [], cfg.conv_fac);
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
fprintf('<strong>%s</strong>: Encode Pos duration = %f\n', mfilename, enc.pos.tvec(end) - enc.pos.tvec(1))
fprintf('<strong>%s</strong>: Recall Pos duration = %f\n', mfilename, rec.pos.tvec(end) - rec.pos.tvec(1))

pos = enc.pos; 
pos.tvec = [enc.pos.tvec, rec.pos.tvec]; 
pos.data = [enc.pos.data, rec.pos.data]; 

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

enc.emg = restrict(emg, enc_t(1), enc_t(2)); 
sleep.emg = restrict(emg, enc_t(2), rec_t(1)); 
rec.emg = restrict(emg, rec_t(1), rec_t(2)); 

enc.acc = restrict(acc, enc_t(1), enc_t(2)); 
sleep.acc = restrict(acc, enc_t(2), rec_t(1)); 
rec.acc = restrict(acc, rec_t(1), rec_t(2)); 
% 
% enc.S = restrict(S, enc_t(1), enc_t(2)); 
% sleep.S = restrict(S, enc_t(2), rec_t(1)); 
% rec.S = restrict(S, rec_t(1), rec_t(2)); 
%% generate the hypnogram

csc_temp = sleep.csc;
csc_temp.data = csc_temp.data(1,:); 
csc_temp.label(2:end) = [];

csc_temp.cfg.hdr(2:end) = [];

if exist('acc', 'var')


hypno = dSub_Sleep_screener(csc_temp, sleep.acc.data(4,:), []); 

else

hypno = dSub_Sleep_screener(csc_temp, sleep.emg, []); 
end
close(221)
%% save it all for output

out = [];
out.meta = meta; 
out.S = S;
out.evts = evts;
out.enc_t = enc_t;
out.rec_t = rec_t; 
% out.csc = csc;
% out.pos = pos;
out.hypno = hypno;
out.history.function{1} = mfilename;
out.history.date{1} = date;
out.history.cfg{1} = cfg; 
end



