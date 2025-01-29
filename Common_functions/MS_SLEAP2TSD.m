function [pos, behav] = MS_SLEAP2TSD(fname,fs, model_in,conv_fac, plot_flag)
%% MS_SLEAP2TSD: loads and collects all SLEAP HD5 files in a directory. 
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
    fname = []; % just use current dir.
    fs = [];
    model_in = [];
    conv_fac = [1 1];
    plot_flag = 0;
elseif nargin == 1
    fs = [];
    model_in = [];
    conv_fac = [1 1];
    plot_flag = 0;
elseif nargin ==2
    model_in = [];
    conv_fac = [1 1];
    plot_flag = 0;
elseif nargin ==3
    conv_fac = [1 1];
    plot_flag = 0;
elseif nargin ==4
    plot_flag = 0;
end


% convert to a cell array if only one file. 
if ~iscell(fname)
    temp = fname; 
    fname = []; 
    fname{1} = temp; 
    clear t; 
end

%%   cycle through all the files and collect the data
data_out = []; 
data_mat = [];
tvec = []; 
vid_ind = []; 
for iF  = 1:length(fname)

    occ_mat = h5read(fname{iF},'/track_occupancy');
    trk_mat = h5read(fname{iF},'/tracks');
    nodes= h5read(fname{iF},'/node_names');
    score_mat = h5read(fname{iF},'/tracking_scores');
    vid_ind{iF} = h5read(fname{iF},'/video_path');  % get the video path.
    
    
    % get the tvec from either the fs argin or extract it from the raw
    % video (slow). 
    if isempty(fs)
        v_obj = VideoReader(parts{end});
        fs = v_obj.FrameRate;
        
        f_c = 0; tvec = [];
        while hasFrame(v_obj)
            readFrame(v_obj);
            f_c = f_c+1;
            tvec{iF}(f_c) = v_obj.CurrentTime;
        end
        clear v_obj;
    else 
        tvec{iF} = 0:length(occ_mat)-1;
        tvec{iF} = tvec{iF}./30;
    end
    
    
    for iN = length(nodes):-1:1
        f_idx = find(contains(nodes, nodes{iN}));
        data_mat{iF, iN} = [trk_mat(:,f_idx,1)./conv_fac(1), trk_mat(:,f_idx,2)./conv_fac(1)]; % grab the data and apply the conversion factor. 
    end
end
   
   
    %make an empty field for each.
    fprintf('Found %i nodes: ', length(nodes));
    for iN = 1:length(nodes)
        nodes{iN}(~isletter(nodes{iN})) = [];  % 
        data_out.(nodes{iN}) = [];
        fprintf('<strong>%s</strong> ', nodes{iN});
    end
    data_out.tvec = []; 
    fprintf('\n');
    
    % put them all together
    fnum = [];
    for iF  = 1:length(fname)
        for iN = 1:length(nodes)
            data_out.(nodes{iN}) = [ data_out.(nodes{iN});  data_mat{iF, iN}];
        end
        data_out.tvec = [data_out.tvec, tvec{iF}]; 
        fnum = [fnum , 1:length(tvec{iF})];
    end

    %% use the real timestamps if they exist. 
    
if exist('timestamps.csv', 'file')
   
        t_file = dir('timestamps.csv');
    
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
        this_cam_idx = nearest_idx(length(data_out.(nodes{1})), TS_len); % get the nearest timestamps.
        data_out.tvec = TS{1}.system_clock{this_cam_idx};
        data_out.tvec(1) = 0;
        frameNum = TS{1}.framenumber{this_cam_idx};
        maxBufferUsed = TS{1}.buffer{this_cam_idx};
        
    end
    data_out.tvec = data_out.tvec./1000; % convert to seconds
    
end

%% apply simple smoothing over any points that are larger than a 3sd jump.
% for iN = 1:length(nodes)
%     nan_idx =abs(zscore(diff(data_out.(nodes{iN})(:,1)))) > 2;
%     
%     data_out.(nodes{iN})(nan_idx,1) = NaN;
%     data_out.(nodes{iN})(:,1) = fillmissing(data_out.(nodes{iN})(:,1), 'spline');
%     nan_idx = abs(zscore(diff(data_out.(nodes{iN})(:,2)))) > 2;
%     data_out.(nodes{iN})(nan_idx,2) = NaN;
%     data_out.(nodes{iN})(:,2) = fillmissing(data_out.(nodes{iN})(:,2), 'spline');
% end

%% compute some other measures.


% get the HD
if length(nodes) ==1 && sum(contains(nodes, 'LED'))>0
    ear_mid(:, 1) = data_out.LED(:,1); 
    ear_mid(:,2) = data_out.LED(:,2); 
    
    HD = nan(length(ear_mid), 1); 

elseif sum(contains(nodes, 'R')) > 0 && sum(contains(nodes, 'L')) > 0 && sum(contains(nodes, 'LED')) > 0
    
    ear_mid(:,1) = (data_out.(nodes{contains(nodes, 'R')})(:,1) + data_out.(nodes{contains(nodes, 'L') & ~contains(nodes, 'LED')})(:,1))/2; % get x mid
    ear_mid(:,2) = (data_out.(nodes{contains(nodes, 'R')})(:,2) + data_out.(nodes{contains(nodes, 'L') & ~contains(nodes, 'LED')})(:,2))/2; % get x mid
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.LED(:,2),ear_mid(:,1) - data_out.LED(:,1)));
    
elseif sum(contains(nodes, 'body', 'IgnoreCase', true))>0 && sum(contains(nodes, 'head', 'IgnoreCase', true))>0
    b = nodes{contains(nodes, 'body', 'IgnoreCase', true)}; 
    h = nodes{contains(nodes, 'head', 'IgnoreCase', true)};
    ear_mid(:,1) = data_out.(b)(:,1);
    ear_mid(:,2) = data_out.(b)(:,2);

    HD = rad2deg(atan2(ear_mid(:,2) - data_out.(h)(:,2),ear_mid(:,1) - data_out.(h)(:,1)));

    
elseif sum(contains(nodes, 'LED'))>0 && sum(contains(nodes, 'Tail'))>0
    ear_mid(:,1) = data_out.Tail(:,1);
    ear_mid(:,2) = data_out.Tail(:,2);
    
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.LED(:,2),ear_mid(:,1) - data_out.LED(:,1)));
    
elseif sum(contains(nodes, 'Body'))>0 && sum(contains(nodes, 'LED'))>0
    ear_mid(:,1) = data_out.Body(:,1);
    ear_mid(:,2) = data_out.Body(:,2);
    
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.LED(:,2),ear_mid(:,1) - data_out.LED(:,1)));
    
    
elseif sum(contains(nodes, 'Green'))>0 && sum(contains(nodes, 'Red'))>0
    ear_mid(:,1) = (data_out.Red(:,1) + data_out.Green(:,1))/2; % get x mid
    ear_mid(:,2) = (data_out.Red(:,2) + data_out.Green(:,2))/2; % get x mid
    
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.Body(:,2),ear_mid(:,1) - data_out.Body(:,1)));
    
elseif sum(contains(nodes, 'nose'))>0 && sum(contains(nodes, 'body'))>0
    ear_mid(:,1) = data_out.body(:,1);
    ear_mid(:,2) = data_out.body(:,2);
    
    HD = rad2deg(atan2(ear_mid(:,2) - data_out.nose(:,2),ear_mid(:,1) - data_out.nose(:,1)));

end

% get the speed from the mid-ear position
vx = dxdt(data_out.tvec,ear_mid(:,1), 'verbose', false);
vy = dxdt(data_out.tvec,ear_mid(:,2), 'verbose', false);




%% convert to pos tsd
data = [];
for iD = 1:length(nodes)
    data = [data, data_out.(nodes{iD})(:,1:2)];
    labels{iD} = nodes{iD};
end

pos = tsd(data_out.tvec,[data, sqrt(vx.^2+vy.^2)', HD]', [labels 'Speed', 'HD']);

if conv_fac(1) == 1 && conv_fac(2) == 1
    pos.units = 'px';
else
    pos.units = 'cm';
end

pos.cfg.conv_fac = conv_fac; 
pos.cfg.fname = fname; 
pos.cfg.vid = vid_ind; 
%% grab a video for a mean frame. 
v_fname  = strsplit(pos.cfg.vid{1}, {'\', '/'}); 
v_fname = v_fname{end}; 

if exist(pos.cfg.vid{1}, 'file')

    vidObj = VideoReader(pos.cfg.vid{1});
    % F=read(vidObj);
    frames = randsample(1:vidObj.NumFrames, 10);
    for ii = length(frames):-1:1
        m =double(read(vidObj,ii));
        F(:,:,ii) = m(:,:,1);
    end
    pos.mean_frame=median(F,3);
    % pos.mean_frame=median(F,4);

elseif exist(v_fname, 'file')
    vidObj = VideoReader(v_fname);
    frames = randsample(1:vidObj.NumFrames, 10);
    for ii = length(frames):-1:1
        m =double(read(vidObj,ii));
        F(:,:,ii) = m(:,:,1);
    end

    pos.mean_frame=median(F,3);

else
    pos.mean_frame = [];

end
    



%% convert to behav format
behav = [];
behav.time = data_out.tvec;
behav.dirName = cd;
behav.numFiles = length(fname);
behav.numFrames = data_out.tvec;
behav.vidNum = fnum;
behav.frameNum = 1:length(data_out.tvec);
behav.maxFramesPerFile = 1000;
behav.height = ceil(max(data_out.(nodes{1})(:,1)));
behav.width =  ceil(max(data_out.(nodes{1})(:,2)));
behav.camNumber = 1;
behav.maxBufferUsed =  NaN;
behav.position = data_out.(nodes{1})(:,1:2);
behav.speed = sqrt(vx.^2+vy.^2)';
behav.HD = HD;
behav.json = []; 
behav.conv_fac = conv_fac; 

%% test out the HD.
if plot_flag
    figure(797)
    c_ord = parula(length(nodes)); 
    s_ord = {'o', 's', 'd', '+', 'x', 'o', 's', 'd', '+', 'x'}; 
    clf
    hold on
    for iF =1:floor(length(pos.data)/20)
        c = 0;
        for ii = 1:length(nodes)
            plot(pos.data(0+ii+c,iF), pos.data(1+ii+c,iF),'.', 'color', c_ord(ii,:), 'marker', s_ord{ii} )
%             plot(pos.data(3,iF), pos.data(4,iF), 'sb')
                    c = c+1;

        end
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
% cd(og_dir);


%% hold for looping version;
%% find all the files

%%%%% TO Do implement model_in catch %%%%%%%%

% 
% og_dir = dir_in;
% cd(dir_in);
% 
% file_list = {};
% d = dir(['*.analysis.h5']);
% rem_idx = zeros(1,length(d));
% for iF = length(d):-1:1
%     parts = strsplit(d(iF).name, {'_', '.'});
%     
%     if isempty(model_in) && any(parts{2} >= '0' & parts{2} <= '9')
%         inter_ver(iF) = str2double(parts{2}(2:end));
%         file_list{iF} = d(iF).name;
%     elseif ~isempty(model_in) && contains(d(iF).name, model_in)
%         inter_ver(iF) = str2double(parts{2}(2:end));
%         file_list{iF} = d(iF).name;
%     else
%         inter_ver(iF) = NaN;
%         file_list{iF} = NaN;
%         rem_idx(iF) = iF;
%     end
% end
% 
% 
% newest_inter_ver = max(inter_ver);
% % loop and find only the DLC versions that use the best trained model (ie most
% % iterations).
% rem_idx = zeros(1,length(file_list));
% for iF = length(file_list):-1:1
%     parts = strsplit(file_list{iF}, {'_', '.'});
%     if str2double(parts{2}(2:end)) ~= newest_inter_ver
%         rem_idx(iF) = iF;
%     end
% %     parts = strsplit(file_list{iF}, {'.'});
% %     
% %     if strcmp(parts{end-2}(1:3), '202'); % check if the output does not sta
% %     Vid_name{iF}(1:strfind(Vid_name{iF},'202'))
%     
% end
% rem_idx(rem_idx == 0) = [];
% file_list(rem_idx) = [];
% 
% %%   cycle through all the files and collect the data
% data_out = [];
% this_node = [];
% 
% for iF  = 1:length(file_list)
%     
%     
%     this_data_out = [];
%     
%     occ_mat = h5read(file_list{iF},'/track_occupancy');
%     trk_mat = h5read(file_list{iF},'/tracks');
%     nodes= h5read(file_list{iF},'/node_names');
%     score_mat = h5read(file_list{iF},'/tracking_scores');
%     vid_ind = h5read(file_list{iF},'/video_path');  % get the video path.
%     
%     parts = strsplit(vid_ind, '/');
%     
%     if isempty(fs)
%         v_obj = VideoReader(parts{end});
%         fs = v_obj.FrameRate;
%         
%         f_c = 0; tvec = [];
%         while hasFrame(v_obj)
%             readFrame(v_obj);
%             f_c = f_c+1;
%             tvec(f_c) = v_obj.CurrentTime;
%         end
%         data_out.tvec = tvec;
%         clear tvec v_obj;
%     else
%         data_out.tvec = 0:length(occ_mat)-1;
%         data_out.tvec = tvec./30;
%     end
%     
%     
%     
%     % get a confidence level for the data
%     
%     
%     
%     %     fields = unique(DLC_labels); % get the parts
%     for iN = 1:length(nodes)
%         f_idx = find(contains(nodes, nodes{iN}));
%         this_node{iF, iN} = [trk_mat(:,f_idx,1), trk_mat(:,f_idx,2)];
%     end
%     
%     
%     
% end
% 
% % apply simple smoothing over any points that are larger than a 3sd jump.
% for iN = 1:length(nodes)
%     nan_idx =abs(zscore(diff(data_out.(nodes{iN})(:,1)))) > 2;
%     
%     data_out.(nodes{iN})(nan_idx,1) = NaN;
%     data_out.(nodes{iN})(:,1) = fillmissing(data_out.(nodes{iN})(:,1), 'spline');
%     nan_idx = abs(zscore(diff(data_out.(nodes{iN})(:,2)))) > 2;
%     data_out.(nodes{iN})(nan_idx,2) = NaN;
%     data_out.(nodes{iN})(:,2) = fillmissing(data_out.(nodes{iN})(:,2), 'spline');
% end
% 
% 
% % convert to cfg units. Default conv fac is 1, so still in pixels.
% for iN = 1:length(nodes)
%     data_out.(nodes{iN})(:,1)  = data_out.(nodes{iN})(:,1)./ conv_fac(1);
%     data_out.(nodes{iN})(:,2)  = data_out.(nodes{iN})(:,2)./ conv_fac(2);
% end
% 
% 
% nan_idx = isnan(data_out.(nodes{1})(:,1));
