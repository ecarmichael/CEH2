function [pos, behav] = MS_DLC2TSD(dir_in, model_in,conv_fac, plot_flag)
%% MS_DLC2TSD: loads and collects all DLC files in a directory. Will skip over files without a number since DLC saves the interation number in the .csv
%
%
%    Inputs:
%    -  dir_in:  [string]  name of the directory to collect.
%
%    - model_in [string]  string to find in the DLC .csv files.  example:
%    '0DLC_resnet_50_II_HD_modelFeb16shuffle1_300000.csv'
%    '0DLC_resnet_50_II_HD_modelFeb16shuffle1_400000.csv',  if model_in =
%    '400000' if will only keep the files containing '400000'.
%
%    - plot_flag: bool do you want to see a movie of the HD / position?
%
%    Outputs:
%    - pos [struct] position data in the tsd format.
%
%    - behav [struct]  data in the 'behav' format for miniscope analysis.
%    [defaults is 0]
%
%
%
%
% EC 2023-01-28  initial version
%% initialize

if nargin == 0
    dir_in = cd; % just use current dir.
    model_in = [];
    conv_fac = [1 1];
    plot_flag = 0;
elseif nargin == 1
    model_in = [];
    conv_fac = [1 1];
    plot_flag = 0;
elseif nargin ==2
    conv_fac = [1 1];
    plot_flag = 0;
elseif nargin ==3
    plot_flag = 0;
end


%% find all the files

%%%%% TO Do implement model_in catch %%%%%%%%


og_dir = dir_in;
cd(dir_in);

file_list = {};
d = dir(['*DLC*.csv']);
rem_idx = zeros(1,length(d));
for iF = length(d):-1:1
    parts = strsplit(d(iF).name, {'_', '.'});
    
    if isempty(model_in) && any(parts{end-1} >= '0' & parts{end-1} <= '9')
        inter_ver(iF) = str2double(parts{end-1});
        file_list{iF} = d(iF).name;
    elseif ~isempty(model_in) && contains(d(iF).name, model_in)
        inter_ver(iF) = str2double(parts{end-1});
        file_list{iF} = d(iF).name;
    else
        inter_ver(iF) = NaN;
        file_list{iF} = NaN;
        rem_idx(iF) = iF;
    end
    f_num(iF) = str2double( regexp(d(iF).name, '\d+','match', 'once' )); % used to sort the files.
end

[~, idx] = sort(f_num);              % sort the file names based on
file_list  = file_list(idx);
inter_ver  = inter_ver(idx);
rem_idx = rem_idx(idx);

% remove empty cells
rem_idx(rem_idx == 0) = [];
file_list(rem_idx) = [];

newest_inter_ver = max(inter_ver);
% loop and find only the DLC versions that use the best trained model (ie most
% iterations).
rem_idx = zeros(1,length(file_list));
for iF = length(file_list):-1:1
    parts = strsplit(file_list{iF}, {'_', '.'});
    if str2double(parts{end-1}) ~= newest_inter_ver
        rem_idx(iF) = iF;
    end
end
rem_idx(rem_idx == 0) = [];
file_list(rem_idx) = [];

%%   cycle through all the files and collect the data

this_field = [];

for iF  = 1:length(file_list)
    DLC_data = table2array(readtable(file_list{iF}, 'Headerlines', 3));
    tbl = readtable(file_list{iF},'Headerlines', 1);
    DLC_labels = tbl.Properties.VariableNames(2:end);
    fields= DLC_labels;
    % hack to have it work for Bryan's data
    if sum(contains(DLC_labels, 'body_center'))>0
        bd_idx = find(contains(DLC_labels, 'body_center')); 
        DLC_labels{bd_idx(1)} = 'body';
                fields{bd_idx(1)} = 'body';

    end
    fields(contains(DLC_labels, '_')) = [];
    %     fields = unique(DLC_labels); % get the parts
    
    for iFields = 1:length(fields)
        f_idx = find(contains(DLC_labels, fields{iFields}));
        this_field{iF, iFields} = [DLC_data(:,f_idx(1)+1), DLC_data(:,f_idx(2)+1), DLC_data(:,f_idx(3)+1)];
    end
end

%make an empty field for each.
fprintf('Found %i fields: ', length(fields));
for iFields = 1:length(fields)
    data_out.(fields{iFields}) = [];
    fprintf('<strong>%s</strong> ', fields{iFields});
end
fprintf('\n');

% put them all together
fnum = [];
for iF  = 1:length(file_list)
    for iFields = 1:length(fields)
        data_out.(fields{iFields}) = [ data_out.(fields{iFields});  this_field{iF, iFields}];
    end
    fnum = [fnum , 1:length(this_field{iF, 1})];
end


% apply simple smoothing over any points that are larger than a 3sd jump.
for iFields = 1:length(fields)
    nan_idx = zscore(abs(diff(data_out.(fields{iFields})(:,1)))) > 3;
    data_out.(fields{iFields})(nan_idx,1) = NaN;
    data_out.(fields{iFields})(:,1) = fillmissing(data_out.(fields{iFields})(:,1), 'spline');
    nan_idx = zscore(abs(diff(data_out.(fields{iFields})(:,2)))) > 3;
    data_out.(fields{iFields})(nan_idx,2) = NaN;
    data_out.(fields{iFields})(:,2) = fillmissing(data_out.(fields{iFields})(:,2), 'spline');
end


% convert to cfg units. Default conv fac is 1, so still in pixels.
for iFields = 1:length(fields)
    data_out.(fields{iFields})(:,1)  = data_out.(fields{iFields})(:,1)./ conv_fac(1);
    data_out.(fields{iFields})(:,2)  = data_out.(fields{iFields})(:,2)./ conv_fac(2);
end


nan_idx = isnan(data_out.(fields{1})(:,1));

%% grab the timestamps.
if ~isempty(dir('*.nvt'))
    fprintf('<strong>%s</strong>: NVT file found, using it for timestamps\n', mfilename)
    
    nvt_f =  dir('*.nvt');
    
    [data_out.tvec, ~, ~, ~, ~, ~, ~] = Nlx2MatVT(nvt_f.name, [1 1 1 1 1 1 ], 1, 1, [] );
    data_out.tvec = data_out.tvec*10^-6; % convert to seconds
    frameNum = 1:length(data_out.tvec); 
    maxBufferUsed = NaN; 
    
else
    
    t_file = dir('time*');
    
    if contains(t_file.name, '.csv')
        TS = readtable(t_file.name);
        data_out.tvec = table2array(TS(:,2));
        frameNum = TS(:,1);
        maxBufferUsed =  max(table2array(TS(:,3)));
        
    elseif contains(t_file.name, '.dat')
        cfg_ts = [];
        cfg_ts.correct_time = 0;
        
        TS = MS_Load_TS(cfg_ts);
        for ii = length(TS{1}.system_clock):-1:1
            TS_len(ii) = length(TS{1}.system_clock{ii});
        end
        this_cam_idx = nearest_idx(length(data_out.(fields{1})), TS_len); % get the nearest timestamps.
        data_out.tvec = TS{1}.system_clock{this_cam_idx};
        data_out.tvec(1) = 0;
        frameNum = TS{1}.framenumber{this_cam_idx};
        maxBufferUsed = TS{1}.buffer{this_cam_idx};
        
    end
    data_out.tvec = data_out.tvec./1000; % convert to seconds
    
    
end
if length(data_out.tvec) ~= length(data_out.(fields{1}))
    fprintf('DLC samples (%.0f) differ from timestamps.dat (%.0f). Trimming...\n',length(data_out.(fields{1})), length(data_out.tvec));
    
    if length(data_out.tvec)< length(data_out.(fields{1}))
        for iFields = 1:length(fields)
            data_out.(fields{iFields}) = data_out.(fields{iFields})(1:length(data_out.tvec),:);
        end
    end
end
%% get the meta data from the json to pull the recording start time.
if ~isempty(dir('*.json'))
    
    j_files = dir('*.json');
    
    fid = fopen([fileparts(j_files.folder) filesep j_files.name]);
    raw = fread(fid,inf);
    str = char(raw');
    Exp_json = jsondecode(str);
    fclose(fid);
    
else
    Exp_json = [];
end

%% compute some other measures.



% get the HD

if length(fields) ==1 && sum(contains(fields, 'LED'))>0
    ear_mid(:, 1) = data_out.LED(:,1); 
    ear_mid(:,2) = data_out.LED(:,2); 
    
    HD = nan(length(ear_mid), 1); 

elseif sum(contains(fields, 'R')) > 0 && sum(contains(fields, 'L')) > 0 && sum(contains(fields, 'LED')) > 0
    
    ear_mid(:,1) = (data_out.(fields{contains(fields, 'R')})(:,1) + data_out.(fields{contains(fields, 'L') & ~contains(fields, 'LED')})(:,1))/2; % get x mid
    ear_mid(:,2) = (data_out.(fields{contains(fields, 'R')})(:,2) + data_out.(fields{contains(fields, 'L') & ~contains(fields, 'LED')})(:,2))/2; % get x mid
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.LED(:,2),ear_mid(:,1) - data_out.LED(:,1)));
    
    
elseif sum(contains(fields, 'body'))>0 && sum(contains(fields, 'head'))>0
    ear_mid(:,1) = data_out.body(:,1);
    ear_mid(:,2) = data_out.body(:,2);
    
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.head(:,2),ear_mid(:,1) - data_out.head(:,1)));
    
    
elseif sum(contains(fields, 'Green'))>0 && sum(contains(fields, 'Red'))>0
    ear_mid(:,1) = (data_out.Red(:,1) + data_out.Green(:,1))/2; % get x mid
    ear_mid(:,2) = (data_out.Red(:,2) + data_out.Green(:,2))/2; % get x mid
    
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.Body(:,2),ear_mid(:,1) - data_out.Body(:,1)));
    
elseif sum(contains(fields, 'nose'))>0 && sum(contains(fields, 'body'))>0
    ear_mid(:,1) = data_out.body(:,1);
    ear_mid(:,2) = data_out.body(:,2);
    
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.nose(:,2),ear_mid(:,1) - data_out.nose(:,1)));

end

% get the speed from the mid-ear position
vx = dxdt(data_out.tvec,ear_mid(:,1));
vy = dxdt(data_out.tvec,ear_mid(:,2));




%% convert to pos tsd
data = [];
for iD = 1:length(fields)
    data = [data, data_out.(fields{iD})(:,1:2)];
    labels{iD} = fields{iD};
end

pos = tsd(data_out.tvec,[data, sqrt(vx.^2+vy.^2)', HD]', [labels 'Speed', 'HD']);

if conv_fac(1) == 1 && conv_fac(2) == 1
    pos.units = 'px';
else
    pos.units = 'cm';
end

pos.cfg.json = Exp_json;

if isfield(pos.cfg.json, 'recordingStartTime') && ~isfield(pos.cfg.json, 'msecSinceEpoch')
    pos.cfg.json.msecSinceEpoch = pos.cfg.json.recordingStartTime.msecSinceEpoch;
end

%% convert to behav format
behav = [];
behav.time = data_out.tvec;
behav.dirName = cd;
behav.numFiles = length(file_list);
behav.numFrames = data_out.tvec;
behav.vidNum = fnum;
behav.frameNum = frameNum;
behav.maxFramesPerFile = 1000;
behav.height = ceil(max(data_out.(fields{1})(:,1)));
behav.width =  ceil(max(data_out.(fields{1})(:,2)));
behav.camNumber = 1;
behav.maxBufferUsed =  maxBufferUsed;
behav.position = data_out.(fields{1})(:,1:2);
behav.speed = sqrt(vx.^2+vy.^2)';
behav.HD = HD;
behav.json = Exp_json;

%% test out the HD.
if plot_flag
    figure(797)
    clf
    hold on
    for iF =1:floor(length(pos.data)/20)
        plot(pos.data(1,iF), pos.data(2,iF), 'or')
        plot(pos.data(3,iF), pos.data(4,iF), 'sb')
        plot([pos.data(3,iF), pos.data(1,iF)],[pos.data(4,iF), pos.data(2,iF)], 'k')
        
        %     plot(behav.position(iF,1), behav.position(iF,2), 'or')
        %     plot([data_out.R_ear(iF,1),  data_out.L_ear(iF,1)], [data_out.R_ear(iF,2),  data_out.L_ear(iF,2)], 'sr')
        %     plot(ear_mid(iF, 1), ear_mid(iF,2), '.b')
        %     plot(data_out.LED(iF,1), data_out.LED(iF,2), 'xg')
        %     plot([ear_mid(iF,1), behav.position(iF,1)],[ear_mid(iF,2), behav.position(iF,2)], 'k')
        xlim([min(pos.data(1,:)) max(pos.data(1,:))]);
        ylim([min(pos.data(2,:)) max(pos.data(2,:))]);
        text(min(pos.data(1,:)), min(pos.data(2,:))+((max(pos.data(2,:)) - min(pos.data(2,:)))/18),['HD (deg): ' num2str(pos.data(end,iF),3)]);
        text(min(pos.data(1,:)), min(pos.data(2,:))+((max(pos.data(2,:)) - min(pos.data(2,:)))/10),['Speed (' pos.units '/s): ' num2str(pos.data(end-1,iF),3)]);
        drawnow
        pause(0.0303/300)
        cla(gca)
    end
end
cd(og_dir);
