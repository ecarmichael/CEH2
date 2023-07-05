function [CSC, S, pos, evt] = Master_NOL_preprocess(cfg_in, Kilo_dir, OE_dir, DLC_dir, save_dir)
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
cfg_def.conv_fac = [6.6  5.8];
cfg_def.behav_label = 'Body';
cfg_def.TTL_trig = '6'; 
cfg_def.TTL_frame = '4'; 
cfg_def.csc_chan = {'CH44', 'CH51'};
cfg_def.emg_chan = {64}; 


cfg = ProcessConfig(cfg_def, cfg_in);


CSC = [];
S = [];
pos = [];
evt = [];


%%  Get the position data if it is there

if ~isempty(DLC_dir)
    [pos, behav] = MS_DLC2TSD(DLC_dir, [], cfg.conv_fac);
end
fprintf('<strong>%s</strong>: Pos duration = %f\n', mfilename, pos.tvec(end) - pos.tvec(1))

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


%% align the behav times to OE times


pause

end



