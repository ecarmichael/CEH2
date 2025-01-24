function [pos, behav] = MS_DLC2TSD_single(fname, vname, conv_fac, plot_flag)
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

if nargin == 1
    error('needs a video file for time');
elseif nargin == 2
    conv_fac = [1 1];
    plot_flag = 0;
elseif nargin ==3
    plot_flag = 0;
end


%% find all the files

% get the name of the corresponding video. 





%%   cycle through all the files and collect the data

this_field = [];

    DLC_data = table2array(readtable(fname, 'Headerlines', 3));
    tbl = readtable(fname,'Headerlines', 1);
    DLC_labels = tbl.Properties.VariableNames(2:end);
    fields= DLC_labels;
    % hack to have it work for Bryan's data
%     if sum(contains(DLC_labels, 'body_center'))>0
%         bd_idx = find(contains(DLC_labels, 'body_center')); 
%         DLC_labels{bd_idx(1)} = 'body';
%                 fields{bd_idx(1)} = 'body';
%     end

% catch extra fields with numbers at the end. 
rm_idx = []; 
for f = length(fields):-1:1
    if ~isnan(str2double(fields{f}(end)))
        rm_idx(f) = true;
    else
        rm_idx(f) = false; 
    end
    
end

fields(logical(rm_idx)) = []; 
    
%     fields(contains(DLC_labels, '_')) = [];
    
    %     fields = unique(DLC_labels); % get the parts
    
    for iFields = 1:length(fields)
        f_idx = find(contains(DLC_labels, fields{iFields}));
        this_field{iFields} = [DLC_data(:,f_idx(1)+1), DLC_data(:,f_idx(2)+1), DLC_data(:,f_idx(3)+1)];
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
    for iFields = 1:length(fields)
        data_out.(fields{iFields}) = [ data_out.(fields{iFields});  this_field{iFields}];
    end
    fnum = [fnum , 1:length(this_field{1})];



% apply simple smoothing over any points that are larger than a 3sd jump.
for iFields = 1:length(fields)
    nan_idx =abs(zscore(diff(data_out.(fields{iFields})(:,1)))) > 2;
    
    data_out.(fields{iFields})(nan_idx,1) = NaN;
    data_out.(fields{iFields})(:,1) = fillmissing(data_out.(fields{iFields})(:,1), 'spline');
    nan_idx = abs(zscore(diff(data_out.(fields{iFields})(:,2)))) > 2;
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
    
    
elseif contains(vname, '.mp4')
    
    v = VideoReader(vname);
    data_out.tvec = (0:v.NumFrames)./v.FrameRate; 
    frameNum = v.NumFrames;
    
    
    
    
else
    
    t_file = dir('time*');
    
    if contains(t_file.name, '.csv')
        TS = readtable(t_file.name);
        data_out.tvec = table2array(TS(:,2));
        frameNum = table2array(TS(:,1));
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
    elseif length(data_out.tvec)> length(data_out.(fields{1}))
        data_out.tvec = data_out.tvec(1:length(data_out.(fields{iFields}))); 
        
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
    
elseif sum(contains(fields, 'Body'))>0 && sum(contains(fields, 'Neck'))>0
    ear_mid(:,1) = data_out.Body(:,1);
    ear_mid(:,2) = data_out.Body(:,2);
    
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.Neck(:,2),ear_mid(:,1) - data_out.Neck(:,1)));
    
elseif sum(contains(fields, 'LED'))>0 && sum(contains(fields, 'Tail'))>0
    ear_mid(:,1) = data_out.Tail(:,1);
    ear_mid(:,2) = data_out.Tail(:,2);
    
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.LED(:,2),ear_mid(:,1) - data_out.LED(:,1)));
    
elseif sum(contains(fields, 'Body'))>0 && sum(contains(fields, 'LED'))>0
    ear_mid(:,1) = data_out.Body(:,1);
    ear_mid(:,2) = data_out.Body(:,2);
    
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.LED(:,2),ear_mid(:,1) - data_out.LED(:,1)));
    
    
elseif sum(contains(fields, 'Green'))>0 && sum(contains(fields, 'Red'))>0
    ear_mid(:,1) = (data_out.Red(:,1) + data_out.Green(:,1))/2; % get x mid
    ear_mid(:,2) = (data_out.Red(:,2) + data_out.Green(:,2))/2; % get x mid
    
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.Body(:,2),ear_mid(:,1) - data_out.Body(:,1)));
    
elseif sum(contains(fields, 'nose'))>0 && sum(contains(fields, 'body'))>0
    ear_mid(:,1) = data_out.body(:,1);
    ear_mid(:,2) = data_out.body(:,2);
    
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.nose(:,2),ear_mid(:,1) - data_out.nose(:,1)));
elseif sum(contains(fields, 'Ear_R'))>0 && sum(contains(fields, 'Ear_L'))>0
       ear_mid(:,1) = (data_out.Ear_R(:,1) + data_out.Ear_L(:,1))/2; % get x mid
        ear_mid(:,2) = (data_out.Ear_R(:,2) + data_out.Ear_L(:,2))/2; % get x mid
        HD = rad2deg(atan2(ear_mid(:,2) - data_out.Body(:,2),ear_mid(:,1) - data_out.Body(:,1)));

    
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
pos.cfg.conv_fac = conv_fac; 

if isfield(pos.cfg.json, 'recordingStartTime') && ~isfield(pos.cfg.json, 'msecSinceEpoch')
    pos.cfg.json.msecSinceEpoch = pos.cfg.json.recordingStartTime.msecSinceEpoch;
end

%% grab a video for a mean frame. 
avi_dir = dir([ vname]);
if ~isempty(avi_dir)
    vidObj = VideoReader([avi_dir(1).folder filesep avi_dir(1).name]);   
    
    F=read(vidObj, [1,200]);                             
    pos.mean_frame=median(F,4);
    clear vidObj F
else
    pos.mean_frame = []; 
end

%% convert to behav format
behav = [];
behav.time = data_out.tvec;
behav.dirName = cd;
behav.numFiles = 1;
behav.numFrames = data_out.tvec;
behav.vidNum = fnum;
behav.frameNum = frameNum;
behav.maxFramesPerFile = 1000;
behav.height = ceil(max(data_out.(fields{1})(:,1)));
behav.width =  ceil(max(data_out.(fields{1})(:,2)));
behav.camNumber = 1;
behav.maxBufferUsed =  NaN;
behav.position = data_out.(fields{1})(:,1:2);
behav.speed = sqrt(vx.^2+vy.^2)';
behav.HD = HD;
behav.json = Exp_json;
behav.conv_fac = conv_fac; 

%% test out the HD.
if plot_flag
    figure(797)
    vidObj = VideoReader(vname); 
    clf
    hold on
    for iF =1:floor(length(pos.data)/20)
        
        imagesc(read(vidObj,iF))
        
        plot(pos.data(1,iF)*conv_fac(1), pos.data(2,iF)*conv_fac(2), 'or')
        plot(pos.data(3,iF)*conv_fac(1), pos.data(4,iF)*conv_fac(2), 'sb')
        plot(ear_mid(iF,1)*conv_fac(1), ear_mid(iF,2)*conv_fac(2), 'dy')
        plot([pos.data(3,iF), pos.data(1,iF)]*conv_fac(1),[pos.data(4,iF), pos.data(2,iF)]*conv_fac(2), 'k')
        
        %     plot(behav.position(iF,1), behav.position(iF,2), 'or')
        %     plot([data_out.R_ear(iF,1),  data_out.L_ear(iF,1)], [data_out.R_ear(iF,2),  data_out.L_ear(iF,2)], 'sr')
        %     plot(ear_mid(iF, 1), ear_mid(iF,2), '.b')
        %     plot(data_out.LED(iF,1), data_out.LED(iF,2), 'xg')
        %     plot([ear_mid(iF,1), behav.position(iF,1)],[ear_mid(iF,2), behav.position(iF,2)], 'k')
        xlim([min(pos.data(1,:)*conv_fac(1)) max(pos.data(1,:)*conv_fac(1))]);
        ylim([min(pos.data(2,:)*conv_fac(2)) max(pos.data(2,:)*conv_fac(2))]);
        text(min(pos.data(1,:)*conv_fac(1)), min(pos.data(2,:))+((max(pos.data(2,:)*conv_fac(2)) - min(pos.data(2,:)*conv_fac(2)))/18),['HD (deg): ' num2str(pos.data(end,iF),3)]);
        text(min(pos.data(1,:)*conv_fac(1)), min(pos.data(2,:))+((max(pos.data(2,:)*conv_fac(2)) - min(pos.data(2,:)*conv_fac(2)))/10),['Speed (' pos.units '/s): ' num2str(pos.data(end-1,iF),3)]);
        drawnow
        pause(0.0303/300)
        cla(gca)
    end
end
% cd(og_dir);
