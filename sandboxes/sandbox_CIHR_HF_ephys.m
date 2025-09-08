%% load some OE data

if ispc
    data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\HF\HF_1_2025-09-05_15-53-34_TFC_REC\Record Node 118';
elseif ismac
    data_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/CIHR_2025/HF/HF_1_2025-09-05_15-53-34_TFC_REC/Record Node 118';
    oe_dir = '/Users/ecar/Documents/Github/open_ephys_matlab_tools';

end
addpath(genpath(oe_dir))
%

evts = OE_LoadEvents();

csc_list = dir('*.continuous');
csc= []; labels = [];
for ii = 1:length(csc_list)

    if ii == 1
        [data, tvec, info] = load_open_ephys_data(csc_list(ii).name);

        csc = tsd(tvec, data);
        labels{ii} = info.header.channel;

    else
        [data, ~, info] = load_open_ephys_data(csc_list(ii).name);
        csc.data =[csc.data, data];
        labels{ii} = info.header.channel;

    end


end

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
