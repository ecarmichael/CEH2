function [pos, behav] = MS_DLC2TSD_2(dir_in, model_in,conv_fac, plot_flag,divide_flag)
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
%    - divide_flag: boolean. 1 For create an independant file instead of
%    adding them to a single file
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
% BC 2024-01-30 V2
%% initialize

if nargin == 0
    dir_in = cd; % just use current dir.
    model_in = [];
    conv_fac = [1 1];
    plot_flag = 0;
    divide_flag=1;
elseif nargin == 1
    model_in = [];
    conv_fac = [1 1];
    plot_flag = 0;
    divide_flag=1;
elseif nargin ==2
    conv_fac = [1 1];
    plot_flag = 0;
elseif nargin ==3
    plot_flag = 0;
    divide_flag=1;
elseif nargin ==4
    divide_flag=1;
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

%So far this-field contains a row per file and a coloum for each point-type (field) tracked in DLC 
% Make an empty field for each one.

fprintf('Found %i fields: ', length(fields));
for iFields = 1:length(fields)
    data_out.(fields{iFields}) = [];
    fprintf('<strong>%s</strong> ', fields{iFields});
end
fprintf('\n');
% Lets split this to either have the option to put them together of gather
% an individual structure form each file 

% put them all together
if divide_flag==0
    fnum = [];
    for iF  = 1:length(file_list)
        for iFields = 1:length(fields)
            data_out.(fields{iFields}) = [ data_out.(fields{iFields});  this_field{iF, iFields}];
        end
        fnum = [fnum , 1:length(this_field{iF, 1})];
    end
%Divide them
%All data is a cell array with that stores a structure for each file being analized
else
    all_data=struct();
    fnum = [];
    %all_data=cell(length(file_list),1);
    for iF =   1:length(file_list)
       for iFields = 1:length(fields)
            data_out.(fields{iFields}) = [this_field{iF, iFields}];
       end
       all_data.("File"+iF)=data_out;
       fnum = [fnum , 1:length(this_field{iF, 1})];
    end
end

% apply simple smoothing over any points that are larger than a 3sd jump.
if divide_flag
    for iFiles= 1:length(fieldnames(all_data))
        for iFields = 1:length(fields)
            nan_idx =abs(zscore(diff(all_data.("File"+iFiles).(fields{iFields})(:,1)))) > 2;
            all_data.("File"+iFiles).(fields{iFields})(nan_idx,1) = NaN;
            all_data.("File"+iFiles).(fields{iFields})(:,1) = fillmissing(all_data.("File"+iFiles).(fields{iFields})(:,1), 'spline');

            nan_idx = abs(zscore(diff(all_data.("File"+iFiles).(fields{iFields})(:,2)))) > 2;
            all_data.("File"+iFiles).(fields{iFields})(nan_idx,2) = NaN;
            all_data.("File"+iFiles).(fields{iFields})(:,2) = fillmissing(all_data.("File"+iFiles).(fields{iFields})(:,2), 'spline');
        end
    end
else
    for iFields = 1:length(fields)
        nan_idx =abs(zscore(diff(data_out.(fields{iFields})(:,1)))) > 2;

        data_out.(fields{iFields})(nan_idx,1) = NaN;
        data_out.(fields{iFields})(:,1) = fillmissing(data_out.(fields{iFields})(:,1), 'spline');
        nan_idx = abs(zscore(diff(data_out.(fields{iFields})(:,2)))) > 2;
        data_out.(fields{iFields})(nan_idx,2) = NaN;
        data_out.(fields{iFields})(:,2) = fillmissing(data_out.(fields{iFields})(:,2), 'spline');
    end
end


% convert to cfg units. Default conv fac is 1, so still in pixels.
if divide_flag
    for iFiles=1:length(fieldnames(all_data))
        for iFields = 1:length(fields)
            all_data.("File"+iFiles).(fields{iFields})(:,1)  = all_data.("File"+iFiles).(fields{iFields})(:,1)./ conv_fac(1);
            all_data.("File"+iFiles).(fields{iFields})(:,2)  = all_data.("File"+iFiles).(fields{iFields})(:,2)./ conv_fac(2);
        end
    end

else
    for iFields = 1:length(fields)
        data_out.(fields{iFields})(:,1)  = data_out.(fields{iFields})(:,1)./ conv_fac(1);
        data_out.(fields{iFields})(:,2)  = data_out.(fields{iFields})(:,2)./ conv_fac(2);
    end
    nan_idx = isnan(data_out.(fields{1})(:,1));
end


%% grab the timestamps.

if ~isempty(dir('*.nvt'))
    fprintf('<strong>%s</strong>: NVT file found, using it for timestamps\n', mfilename)
    nvt_f =  dir('*.nvt');
    if divide_flag
        [all_data.tvec, ~, ~, ~, ~, ~, ~] = Nlx2MatVT(nvt_f.name, [1 1 1 1 1 1 ], 1, 1, [] );
        all_data.tvec = all_data.tvec*10^-6; % convert to seconds
        frameNum = 1:length(all_data.tvec);
        maxBufferUsed = NaN;
    else

        [data_out.tvec, ~, ~, ~, ~, ~, ~] = Nlx2MatVT(nvt_f.name, [1 1 1 1 1 1 ], 1, 1, [] );
        data_out.tvec = data_out.tvec*10^-6; % convert to seconds
        frameNum = 1:length(data_out.tvec);
        maxBufferUsed = NaN;
    end
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
    if divide_flag
        all_data.tvec=all_data.tvec./1000;
    else
        data_out.tvec = data_out.tvec./1000; % convert to seconds
    end

end

%Here we want to avoid trimming timestamps
if divide_flag
    temp_evts=LoadEvents([]);
    ievts= length(temp_evts.t{1,1});
    idlc=length(dir(['*DLC*']));
    % Make sure that the number of events correspond to the #videos+1
    %In a normal recording session for Bryn's experiments. You would end with
    %two videos and three events after deleting the intermediate video.
    if idlc==2 && ievts==3;
        fprintf('<strong>The number of evts corresponds to the number of videos for Bryans experiments</strong> \n Using interval 1 and 2 for extracting the timestamps \n');
        idx_start=find(contains(temp_evts.label,'Starting Recording'));
        idx_end=find(contains(temp_evts.label,'Stopping Recording'));
        range=[temp_evts.t{idx_start}(1),temp_evts.t{idx_end}(1);temp_evts.t{idx_start}(3),temp_evts.t{idx_end}(3)];
    else
        fprintf('<strong>The number of evts DO NOT corresponds to the number of videos for Bryans experiments</strong> \n');
        valid= false;
        while ~valid
            prompt=sprintf("Please provide the index of the recording that you want to use for your <strong>first</strong> interval. There are %d intervals in this experiment \n", ievts);
            idx1_str=input(prompt,"s");
            % Check if the input can be converted to a number
            if ~isempty(idx1_str) && all(isstrprop(idx1_str, 'digit'))
                idx1 = str2double(idx1_str); % Convert the input to a number
                % Check if the conversion was successful
                if ~isnan(idx1) && isreal(idx1) && idx1 > 0
                    valid = true; % Set the flag to true to exit the loop
                end
            end
            % If the input is not valid, display a message and ask again
            if ~valid
                disp('Invalid input. Please enter a valid positive integer.');
            end
        end

        valid= false;
        while ~valid
            prompt=sprintf("Please provide the index of the recording that you want to use for your <strong>second</strong> interval. There are %d intervals in this experiment \n", ievts);
            idx2_str=input(prompt,"s");

            % Check if the input can be converted to a number
            if ~isempty(idx2_str) && all(isstrprop(idx2_str, 'digit'))
                idx2 = str2double(idx2_str); % Convert the input to a number
                % Check if the conversion was successful
                if ~isnan(idx2) && isreal(idx2) && idx2 > 0 && idx1_str~=idx2_str
                    valid = true; % Set the flag to true to exit the loop
                end
            end
            % If the input is not valid, display a message and ask again
            if ~valid
                disp('Invalid input. Please enter a valid positive integer that was not your previous input.');
            end
        end
        idx_start=find(contains(temp_evts.label,'Starting Recording'));
        idx_end=find(contains(temp_evts.label,'Stopping Recording'));

        fprintf('Great! This works. The intervals <strong>%d</strong> and <strong>%d</strong> will be used to extract the times\n', idx1,idx2);
        range=[temp_evts.t{idx_start}(idx1),temp_evts.t{idx_end}(idx1);temp_evts.t{idx_start}(idx2),temp_evts.t{idx_end}(idx2)];
        %If the values do macth, then use the values between 'Starting Recording"
        %and 'Stopping Recording' 1 and 3 as reference for creating a tvec
    end
else
    %%% ----To do-----.The following code assumes that the tvec would
    %%% correspond to the number of frames. This is not the case if you use
    if length(data_out.tvec) ~= length(data_out.(fields{1}))
        fprintf('DLC samples (%.0f) differ from timestamps.dat (%.0f). Trimming position data...\n',length(data_out.(fields{1})), length(data_out.tvec));
        %%% This only trims the position data if the tvec vector is smaller than the position-data  vector
        if length(data_out.tvec)< length(data_out.(fields{1}))
            fprintf('Trimming position data...\n');
            for iFields = 1:length(fields)
                data_out.(fields{iFields}) = data_out.(fields{iFields})(1:length(data_out.tvec),:);
            end
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
if divide_flag

    %HD=struct([]);
    for iF= 1:length(file_list)
        ear_mid=[];
        if length(fields) ==1 && sum(contains(fields, 'LED'))>0
            ear_mid(:, 1) = all_data.("File"+iF).LED(:,1);
            ear_mid(:,2) = all_data.("File"+iF).LED(:,2);
            all_data.("File"+iF).HD= nan(length(ear_mid), 1);
            all_data.("File"+iF).vx= dxdt(all_data.tvec,ear_mid(:,1));
            all_data.("File"+iF).vy = dxdt(all_data.tvec,ear_mid(:,2));

        elseif sum(contains(fields, 'R')) > 0 && sum(contains(fields, 'L')) > 0 && sum(contains(fields, 'LED')) > 0
            % Extract the fields
            R_index = contains(fields, 'R');
            L_index = contains(fields, 'L') & ~contains(fields, 'LED'); % Exclude 'LED' containing fields
            LED_index = contains(fields, 'LED');

            % Check if the conditions return valid indices
            if any(R_index) && any(L_index) && any(LED_index)
                % Calculate the x mid
                ear_mid(:,1) = (all_data.("File"+iF).(fields{R_index})(:,1) + all_data.("File"+iF).(fields{L_index})(:,1))/2;
                ear_mid(:,2) = (all_data.("File"+iF).(fields{R_index})(:,2) + all_data.("File"+iF).(fields{L_index})(:,2))/2;

                % Check if the dimensions of ear_mid and LED match
                if isfield(all_data.("File"+iF), 'LED') && isequal(size(all_data.("File"+iF).LED), size(ear_mid))
                    % Compute HD
                    all_data.("File"+iF).HD = rad2deg(atan2(ear_mid(:,2) - all_data.("File"+iF).LED(:,2), ear_mid(:,1) - all_data.("File"+iF).LED(:,1)));
                else
                    disp("LED data is missing or does not match the size of ear_mid.");
                end
            end
            all_data.("File"+iF).vx= dxdt(all_data.tvec,ear_mid(:,1));
            all_data.("File"+iF).vy = dxdt(all_data.tvec,ear_mid(:,2));

        elseif sum(contains(fields, 'body'))>0 && sum(contains(fields, 'head'))>0
            ear_mid(:,1) = all_data.("File"+iF).body(:,1);
            ear_mid(:,2) = all_data.("File"+iF)(:,2);

            all_data.("File"+iF).HD = rad2deg(atan2(ear_mid(:,2) - all_data.("File"+iF).head(:,2),ear_mid(:,1) - all_data.("File"+iF).head(:,1)));
            all_data.("File"+iF).vx= dxdt(all_data.tvec,ear_mid(:,1));
            all_data.("File"+iF).vy = dxdt(all_data.tvec,ear_mid(:,2));


        elseif sum(contains(fields, 'Green'))>0 && sum(contains(fields, 'Red'))>0
            ear_mid(:,1) = (all_data.("File"+iF).Red(:,1) + all_data.("File"+iF).Green(:,1))/2; % get x mid
            ear_mid(:,2) = (all_data.("File"+iF).Red(:,2) + all_data.("File"+iF).Green(:,2))/2; % get x mid

            all_data.("File"+iF).HD = rad2deg(atan2(ear_mid(:,2) - all_data.("File"+iF).Body(:,2),ear_mid(:,1) - all_data.("File"+iF).Body(:,1)));
            all_data.("File"+iF).vx= dxdt(all_data.tvec,ear_mid(:,1));
            all_data.("File"+iF).vy = dxdt(all_data.tvec,ear_mid(:,2));

        elseif sum(contains(fields, 'nose'))>0 && sum(contains(fields, 'body'))>0
            ear_mid(:,1) = all_data.("File"+iF).body(:,1);
            ear_mid(:,2) = all_data.("File"+iF).body(:,2);

            all_data.("File"+iF).HD = rad2deg(atan2(ear_mid(:,2) - all_data.("File"+iF).nose(:,2),ear_mid(:,1) - all_data.("File"+iF).nose(:,1)));
            all_data.("File"+iF).vx= dxdt(all_data.tvec,ear_mid(:,1));
            all_data.("File"+iF).vy = dxdt(all_data.tvec,ear_mid(:,2));

        end
    end

else
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
    vx = dxdt(data_out.tvec,ear_mid(:,1));
    vy = dxdt(data_out.tvec,ear_mid(:,2));
end

%% convert to pos tsd

if divide_flag
    data = struct();
    pos= struct();

    for iF=1:length(file_list)
        temp_data=[];
        temp_pos=[];

        for iD = 1:length(fields)
            temp_data= [temp_data,all_data.("File"+iF).(fields{iD})(:,1:2)];
            labels{iD} = fields{iD};
        end
        data.("File"+iF)= temp_data;
        %Verify the Tvec below, you may need to trim in the steps above
        vx=all_data.("File"+iF).vx;
        vy=all_data.("File"+iF).vy;
        temp_pos = tsd(all_data.tvec,[temp_data, sqrt(vx.^2+vy.^2)', all_data.("File"+iF).HD]', [labels 'Speed', 'HD']);
        pos.("File"+iF) = temp_pos;
    end
    % If you dont want to divide
else

    data = [];
    for iD = 1:length(fields)
        data = [data, data_out.(fields{iD})(:,1:2)];
        labels{iD} = fields{iD};
    end
    pos = tsd(data_out.tvec,[data, sqrt(vx.^2+vy.^2)', HD]', [labels 'Speed', 'HD']);
end

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

%%-----To do----- Implement the divide version of this section of code
% ----------------------------------You are here in the fx
if divide_flag
    behav = [];
    for iF=1:length(file_list)
        behav.time = all_data.tvec;
        behav.("File" +iF).dirName = cd;
        behav.("File" +iF).numFiles = length(file_list);
        behav.("File" +iF).numFrames = all_data.tvec;
        behav.("File" +iF).vidNum = fnum;
        behav.("File" +iF).frameNum = frameNum;
        behav.("File" +iF).maxFramesPerFile = 1000;
        behav.("File" +iF).height = ceil(max(all_data.("File"+iF).(fields{1})(:,1)));
        behav.("File" +iF).width =  ceil(max(all_data.("File"+iF).(fields{1})(:,2)));
        behav.("File" +iF).camNumber = 1;
        behav.("File" +iF).maxBufferUsed =  maxBufferUsed;
        behav.("File" +iF).position = all_data.("File"+iF).(fields{1})(:,1:2);
        behav.("File" +iF).speed = sqrt(all_data.("File"+iF).vx.^2+all_data.("File"+iF).vy.^2)';
        behav.("File" +iF).HD = all_data.("File"+iF).HD;
        behav.("File" +iF).json = Exp_json;
    end
else
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
end
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
