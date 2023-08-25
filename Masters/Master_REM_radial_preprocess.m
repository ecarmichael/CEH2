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



cfg = ProcessConfig(cfg_def, cfg_in);

Encode = [];
Sleep = [];
Recall = [];


%% check for meta
cd(data_dir)

if ~isempty(dir('*meta.m'))
    meta = MS_Load_meta;
else
    MS_Write_meta_Rad(cd);
    meta = MS_Load_meta;

end


%% load the event
cd(data_dir)


evts = LoadEvents([]);

% get the block times
st_idx = find(contains(evts.label, 'Starting Recording')); 
end_idx = find(contains(evts.label, 'Stopping Recording'));
dur = []; keep = zeros(length(evts.t{st_idx})); 

for ii = length(evts.t{st_idx}):-1:1
    if ((evts.t{end_idx}(ii)  - evts.t{st_idx}(ii))/60) < 5 % skip recordings less than 5min. 
        continue
    end
   dur(end+1) =  evts.t{end_idx}(ii)  - evts.t{st_idx}(ii);
   keep(ii) = 1;
   fprintf('Rec %.0f = %.0fmins \n', ii, dur(end)/60)
    
end

if length(dur) > 3
    error('Too many recordings')
end
% which recordings are real. 
rec_idx = find(keep); 

encode_t = [evts.t{st_idx}(rec_idx(1)), evts.t{end_idx}(rec_idx(1))]; 
sleep_t = [evts.t{st_idx}(rec_idx(2)) , evts.t{end_idx}(rec_idx(2))]; % why was this '- 14400' here: [evts.t{end_idx}(2) - 14400 , evts.t{end_idx}(2)]??
recall_t = [evts.t{st_idx}(rec_idx(3)), evts.t{end_idx}(rec_idx(3))]; 

fprintf('<strong>%s</strong>: Encode: %0.2fmin\n', mfilename,(encode_t(2) - encode_t(1))/60);
fprintf('<strong>%s</strong>: Sleep: %0.2fmin\n', mfilename,(sleep_t(2) -sleep_t(1))/60);
fprintf('<strong>%s</strong>: Recall: %0.2fmin\n', mfilename, (recall_t(2) - recall_t(1))/60);

Encode.evts = restrict(evts, encode_t(1), encode_t(2)); 
Sleep.evts = restrict(evts, sleep_t(1), sleep_t(2));
Recall.evts = restrict(evts, recall_t(1), recall_t(2)); 



%% load the data and trim it
cfg = [];
cfg.getTTnumbers = 0; 
cfg.uint = '64'; %essential for the current NLX offset issue. 
S = LoadSpikes(cfg);


% for ii = length(S.t):-1:1
%     S.t{ii} = S.t{ii}(2:2:end); 
%     
% end

% trim to phases
Encode.S = restrict(S, encode_t(1), encode_t(2)); 
Sleep.S = restrict(S, sleep_t(1), sleep_t(2));
Recall.S = restrict(S, recall_t(1), recall_t(2)); 


%% grab the LFP

cfg = [];
cfg.fc = {meta.goodCSC, meta.goodCSC2, meta.EMG};
cfg.desired_sampling_frequency = 2000;

csc = MS_LoadCSC(cfg);

csc.label{3} = 'emg';

Encode.csc = restrict(csc, encode_t(1), encode_t(2)); 
Sleep.csc = restrict(csc, sleep_t(1), sleep_t(2));
Recall.csc = restrict(csc, recall_t(1), recall_t(2)); 

%% get the position

cfg = [];
% cfg.convFact = [6.4 6.4];

pos = MS_LoadPos(cfg); 


Encode.pos = restrict(pos, encode_t(1), encode_t(2)); 
Sleep.pos = restrict(pos, sleep_t(1), sleep_t(2));
Recall.pos = restrict(pos, recall_t(1), recall_t(2)); 


% correct for pixel 2 cm conversion. 
Encode.pos.data(1,:) = Encode.pos.data(1,:)./6; 
Encode.pos.data(2,:) = Encode.pos.data(2,:)./6; 

Sleep.pos.data(1,:) = Sleep.pos.data(1,:)./4; 
Sleep.pos.data(2,:) = Sleep.pos.data(2,:)./4; 

Recall.pos.data(1,:) = Recall.pos.data(1,:)./6; 
Recall.pos.data(2,:) = Recall.pos.data(2,:)./6; 


% apply light smoothing 
Encode.pos.data(1,:) = smoothdata(Encode.pos.data(1,:),'gaussian', round(1/mode(diff(Encode.pos.tvec)))/2); 
Encode.pos.data(2,:) = smoothdata(Encode.pos.data(2,:),'gaussian', round(1/mode(diff(Encode.pos.tvec)))/2); 

Sleep.pos.data(1,:) = smoothdata(Sleep.pos.data(1,:),'gaussian', round(1/mode(diff(Sleep.pos.tvec)))/2); 
Sleep.pos.data(2,:) = smoothdata(Sleep.pos.data(2,:),'gaussian', round(1/mode(diff(Sleep.pos.tvec)))/2); 

Recall.pos.data(1,:) = smoothdata(Recall.pos.data(1,:),'gaussian', round(1/mode(diff(Recall.pos.tvec)))/2); 
Recall.pos.data(2,:) = smoothdata(Recall.pos.data(2,:),'gaussian', round(1/mode(diff(Recall.pos.tvec)))/2); 

Encode.speed = getLinSpd([],Encode.pos); % linear speed
Sleep.speed = getLinSpd([],Sleep.pos); % linear speed
Recall.speed = getLinSpd([],Recall.pos); % linear speed

% Threshold speed
cfg = []; cfg.method = 'raw'; cfg.operation = 'range'; cfg.threshold = [4 20]; % speed limit in cm/sec
Encode.move.idx = TSDtoIV(cfg,Encode.speed ); % only keep intervals with speed above thresh
Encode.move.cfg = cfg; 

Sleep.move.idx = TSDtoIV(cfg,Sleep.speed ); % only keep intervals with speed above thresh
Sleep.move.cfg = cfg; 

Recall.move.idx = TSDtoIV(cfg,Recall.speed ); % only keep intervals with speed above thresh
Recall.move.cfg = cfg; 


%% get the trials
t_idx = find(ismember(evts.label, 't'));
te_idx = find(ismember(evts.label, 'te'));


t_start =  evts.t{t_idx}; 
t_end = evts.t{te_idx} ; 


% fix overlapping t and te
over_lap_idx = ismember(t_start, t_end);
t_start(over_lap_idx) = [];

t_start_encode = [Encode.csc.tvec(1) t_start(t_start < sleep_t(1))]+60;
t_start_recall = [Recall.csc.tvec(1) t_start(t_start > sleep_t(end))]+60;

t_end_encode = t_end(t_end < sleep_t(1)) ;
t_end_recall = t_end(t_end > sleep_t(end)) ;

Encode.trials = [t_start_encode; t_end_encode];
Encode.ITI = [t_start_encode-60; t_start_encode];
if strcmpi(meta.session, 'Rad3') && strcmpi(meta.subject, 'M30')
    Recall.trials = [t_start_recall(1:4); [t_end_recall Recall.csc.tvec(end)]];
    Recall.ITI = [t_start_recall(1:4)-60; t_start_recall(1:4)];
else
    Recall.trials = [t_start_recall; t_end_recall];
    Recall.ITI = [t_start_recall-60; t_start_recall];
end

% %% generate the hypnogram
% 
% csc_temp = Sleep.csc;
% csc_temp.data = Sleep.csc.data(1,:); 
% csc_temp.label(2:end) = [];
% csc_temp.cfg.hdr(2:end) = [];
% 
% emg_idx = find(contains(lower(Sleep.csc.label), 'emg')); 
% emg_temp = Sleep.csc;
% emg_temp.data = Sleep.csc.data(emg_idx,:); 
% emg_temp.label= [];
% emg_temp.label = Sleep.csc.label(emg_idx); 
% 
% emg_temp.cfg.hdr = [];
% emg_temp.cfg.hdr = Sleep.csc.cfg.hdr(emg_idx); 
% 
% 
% hypno = dSub_Sleep_screener(csc_temp, emg_temp, []); 


%% save it all for output

out = [];
out.meta = meta; 
out.Encode = Encode; 
out.Sleep = Sleep; 
out.Recall = Recall; 

out.history.function{1} = mfilename;
out.history.date{1} = date;
out.history.cfg{1} = cfg; 
end



