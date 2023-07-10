function out = Master_REM_radial_preprocess(cfg_in, data_dir)
%% Master_REM_radial_preprocess:  Extracts the position, spike, and CSC data for a REM radial session.
%
%
%
%    Inputs:
%    - cfg_in: [struct]  configuration parameters to overwrite the defaults
%
%    - dat_dir: [path]  directory with the NLX data
%    Outputs:
%    -
%
%
%
%
% EC 2023-07-10   initial version
%
%
%
%% initialize

cfg_def = [];
cfg_def.conv_fac = [6.6  5.8];
cfg_def.behav_label = 'Body';
cfg_def.TTL_trig = '6'; 
cfg_def.TTL_frame = '4'; 
cfg_def.csc_chan = {'CH44'};
cfg_def.emg_chan = {'CH18'}; 


cfg = ProcessConfig(cfg_def, cfg_in);


CSC = [];
S = [];
pos = [];
evt = [];


%% load the event
cd(OE_dir)

evts = OE_LoadEvents; 

if sum(contains(evts.label, cfg.TTL_trig)) == 1
    evts.label{(contains(evts.label, cfg.TTL_frame))} = 'TTL_frame'; 
    evts.label{(contains(evts.label, cfg.TTL_trig))} = 'TTL_trig'; 
    
    for ii = 1:length(evts.label)
        dur(ii) = evts.t{ii}(end) - evts.t{ii}(1); 
    end
    
    [~, max_TTL] = max(dur); 
    start_t = evts.t{max_TTL}(1);
    end_t = evts.t{max_TTL}(end); 
    
else
    error('No trigger TTL found. Required for synchronizing data')
end

fprintf('<strong>%s</strong>: OE frame TTL (<strong>''%s''</strong>), duration = %.0fsec (%2.1fhr)\n', mfilename,evts.label{max_TTL}, end_t - start_t, (end_t - start_t)/60/60)


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
cfg_csc.fc = cfg.csc_chan;
cfg_csc.desired_sampling_frequency = 2000;
csc = OE_old_csc2TSD(cfg_csc);


cfg_csc = [];
cfg_csc.fc = cfg.emg_chan;
cfg_csc.desired_sampling_frequency = 2000;
emg = OE_old_csc2TSD(cfg_csc);

% 
csc.data(end+1,:) = emg.data;
csc.label{end+1} = 'EMG';
csc.cfg.hdr{end+1} = emg.cfg.hdr{1}; 
%%  Get the position data if it is there

if ~isempty(DLC_dir)
    [pos, ~] = MS_DLC2TSD(DLC_dir, [], cfg.conv_fac);
    if pos.tvec(1) <0
        pos.tvec = pos.tvec + abs(pos.tvec(1)); 
    elseif pos.tvec(1) >0
        pos.tvec  = pos.tvec - abs(pos.tvec(1)); 
    end
end
fprintf('<strong>%s</strong>: Pos duration = %f\n', mfilename, pos.tvec(end) - pos.tvec(1))

%% align the behav times to OE times

fprintf('Pos dur(%.2fsec) - TTLs (%0.2fsec), diff: %0.3f\n', pos.tvec(end) - pos.tvec(1), end_t - start_t, (pos.tvec(end) - pos.tvec(1)) -  (end_t - start_t)) ; 

if ((pos.tvec(end) - pos.tvec(1)) -  (end_t - start_t)) >1
    error('position and TTLs differ by more than 1sec')
end

pos.tvec = pos.tvec + start_t; 


%% generate the hypnogram

csc_temp = csc;
csc_temp.data = csc.data(1,:); 
csc_temp.label(2:end) = [];

csc_temp.cfg.hdr(2:end) = [];


hypno = dSub_Sleep_screener(csc_temp, emg, []); 


%% save it all for output

out = [];
out.S = S;
out.evts = evts;
out.csc = csc;
out.pos = pos;
out.hypno = hypno;
out.history.function{1} = mfilename;
out.history.date{1} = date;
out.history.cfg{1} = cfg; 
end



