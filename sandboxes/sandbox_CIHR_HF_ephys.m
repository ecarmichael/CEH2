%% load some OE data

if ispc
    data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\HF\HF_1_2025-09-08_12-23-48_TFC_REC3_swr\Record Node 119';
    evts_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\HF\HF_1_2025-09-08_12-23-48_TFC_REC3_swr\Record Node 118';
    spk_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\HF\HF_1_2025-09-08_12-23-48_TFC_REC3_swr\Record Node 113\experiment1\recording1\continuous\Intan_RHD_USB-100.Rhythm Data';

elseif ismac
    data_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/CIHR_2025/HF/HF_1_2025-09-05_15-53-34_TFC_REC/Record Node 118';
    oe_dir = '/Users/ecar/Documents/Github/open_ephys_matlab_tools';
addpath(genpath(oe_dir))

end
%
cd(evts_dir)
evts = OE_LoadEvents();


csc_list = dir([data_dir filesep '*.continuous']);
csc= []; labels = [];
for ii = 1:length(csc_list)

    if ii == 1
        [data, tvec, info] = load_open_ephys_data([csc_list(ii).folder filesep csc_list(ii).name]);
        csc = tsd(tvec, data);
        labels{ii} = info.header.channel;
        csc.cfg.hdr{ii}.SamplingFrequency = 30000; 
    else
        [data, ~, info] = load_open_ephys_data([csc_list(ii).folder filesep csc_list(ii).name]);
        csc.data =[csc.data, data];
        labels{ii} = info.header.channel;
        csc.cfg.hdr{ii}.SamplingFrequency = 30000; 
    end
end

csc.data = csc.data'; 
csc.label = labels; 

cfg_in.decimateFactor = 15; 
csc_r = decimate_tsd(cfg_in, csc);

% correct for start of csc. Why is this offset? 
csc_r.tvec = csc_r.tvec- csc_r.tvec(1); 

% loop over events and remove the offset. 
for ii = 1:length(evts.t)
    evts.t{ii} = evts.t{ii} - csc.tvec(1); 
end

%% compare events in binary and OE format. 
cd(evts_dir)
OE_evts = OE_LoadEvents();

cd('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\HF\HF_1_2025-09-08_12-23-48_TFC_REC3_swr\Record Node 113\experiment1\recording1\events\Intan_RHD_USB-100.Rhythm Data\TTL')
BIN_evts = OE_load_binary_evts(cd);

% csc_r = restrict(csc_r, rec_iv);

csc_r.tvec = csc_r.tvec - rec_iv.tstart; 


for ii = 1:size(BIN_evts.t{1},2)
disp(['Binary: ' num2str(BIN_evts.t{1}(1, ii)) ' : ' num2str(BIN_evts.t{1}(2, ii))  '  diff: ' num2str(BIN_evts.t{1}(2, ii)- BIN_evts.t{1}(1, ii)) 's ...' ...
    '(' num2str((BIN_evts.t{1}(2, ii)- BIN_evts.t{1}(1, ii))/60)  'min)'  ])

disp(['OE    : ' num2str(OE_evts.t{1}(1, ii)) ' : ' num2str(OE_evts.t{1}(2, ii))  '  diff: ' num2str(OE_evts.t{1}(2, ii)- OE_evts.t{1}(1, ii)) 's ...' ...
    '(' num2str((OE_evts.t{1}(2, ii)- OE_evts.t{1}(1, ii))/60)  'min)'  ])
end
%% 


rec_iv = iv(evts.t{1}(1), evts.t{1}(2)); 


swr_iv = iv(evts.t{end-1}-.05,evts.t{end-1}+.1); 

swr_iv = MergeIV([], swr_iv); 

swr_r_iv = swr_iv; 
% swr_r_iv = restrict(swr_iv, rec_iv);

for ii = 1:length(swr_r_iv.tstart)
    swr_r_iv.tstart(ii) = swr_r_iv.tstart(ii) - rec_iv.tstart; 
    swr_r_iv.tend(ii) = swr_r_iv.tend(ii) - rec_iv.tstart; 
end

cd(spk_dir)
S = OE_phy2TS(spk_dir);

S_r = S;


for ii = 1:length(S_r.t)
    S_r.t{ii} = S_r.t{ii}-rec_iv.tstart; 
end


ca_idx = []; 
for ii = 1:length(S_r.t)
    p = strsplit(S_r.label{ii}, '-'); 

    if str2double(p{1}) > 75
        ca_idx(ii) = 1; 
    else
        ca_idx(ii) = 0; 
    end
end



%% split CSC and swr_iv in to trials


%% load and make a raster of the spikes

csc_p = csc_r;
csc_p.data = csc_p.data(2,:); 
csc_p.cfg.hdr = []; 
csc_p.cfg.hdr{1} = csc.cfg.hdr{2}; 
% 
% S_r_9 = restrict(S_r, 1970, 2000); 
% csc_p = restrict(csc_p, 1970, 2000); 


% cfg.lfp = csc_p; 

figure(2)
clf

% ax(1) = subplot(4,1,1);
% cfg_plt.display = 'tsd';
% PlotTSDfromIV(cfg_plt, swr_r_iv, csc_p)
% 
% vline(evts.t{2}, 'r')
% vline(evts.t{3}, 'r')

ax(1) = subplot(7,1,1:3);
cfg_mr.lfp = csc_p; 
cfg_mr.evt = swr_r_iv; 
cfg_mr.spkColor = [repmat([0 0 0], 41,1); repmat([0.3467 0.5360  0.6907], 18,1)]; 
cfg_mr.openNewFig = 0; 
h = MultiRaster(cfg_mr, S_r);

ax(2) = subplot(7,1,4:6);
hold on
for ii = 1:size(rate_tsd.data,1)
    plot(rate_tsd.tvec, normalize(rate_tsd.data(ii,:), 'range', [0 1])+ii, 'color', cfg_mr.spkColor(ii,:));

end

% 42-end are shank 4. 

ax(4) = subplot(7,1,7);
cla; 
% cfg_mua = []; 
% cfg_mua.tvec = csc_r.tvec;
% MUA = getMUA(cfg_mua, S_r);
% plot(MUA.tvec, zscore(MUA.data))
% xlim([MUA.tvec(1) MUA.tvec(end)])
plot(rate_tsd.tvec, R(:,keep_idx(1:5)))



vline(evts.t{1}, 'g')
vline(evts.t{2}, 'r')
vline(evts.t{3}, 'b')


linkaxes(ax, 'x')


%% grab a binary channel

chan_idx = 48; 

show = true;
if show
    % Define visualization figure
    set(0,'units','pixels');
    s = get(0,'screensize');
    SCREEN_X = s(3);
    SCREEN_Y = s(4);
    FIGURE_X_SIZE = SCREEN_X / 2;
    FIGURE_Y_SIZE = SCREEN_Y / 3;

    f = figure();
    f.set('Position', [SCREEN_X / 2 0 FIGURE_X_SIZE FIGURE_Y_SIZE]);
end

% Define path to the recording
rec_path = join([DATA_PATH, latest_recording.name]);

% Create a session (loads all data from the most recent recording)
session = Session('/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/CIHR_2025/HF/HF_1_2025-09-05_15-53-34_TFC_REC');

% Get the number of record nodes for this session
nRecordNodes = length(session.recordNodes);

node = session.recordNodes{2};


% 1. Get the first recording
recording = node.recordings{1,1};

% 2. Iterate over all continuous data streams in the recording
streamNames = recording.continuous.keys();


streamName = streamNames{1};
disp(streamName)

% Get the continuous data from the current stream
data = recording.continuous(streamName);

% Plot first channel of continuous data
if show
    plot(data.timestamps, data.samples(chan_idx,:), 'LineWidth', 1.5);
    title(recording.format, recording.format); hold on;
end

% 3. Overlay all available event data
eventProcessors = recording.ttlEvents.keys();
for p = 1:length(eventProcessors)
    processor = eventProcessors{p};
    events = recording.ttlEvents(processor);
    if ~isempty(events)
        for n = 1:length(events.line)
            if show && events.state(n) == 1
                line([events.timestamp(n), events.timestamp(n)], [-4000,2000], 'Color', 'red', 'LineWidth', 0.2);
            end
        end
    end

end
