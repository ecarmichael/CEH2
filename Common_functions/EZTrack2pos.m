function [pos, spd, scale_fac, ROI] =  EZTrack2pos(data_dir, plot_flag)
%% EZTrack2pos: converts the output csv files from ezTrack (https://github.com/denisecailab/ezTrack) to the pos TSD format.
%
%
%
%    Inputs:
%    - data_dir [optional]: string  Path to data
%
%    - plot_flag [optional; default on]: binary - toggles plots
%
%    Outputs:
%    -
%
%
%
%
% EC 2023-01-18   initial version
%
%
%
%% initialize

if nargin < 1
    data_dir = cd;
    plot_flag = 1;
elseif nargin <2
    plot_flag = 1;
end

if isempty(data_dir)
    data_dir = cd;
end

%% find and collect all position files

d = dir('*Output.csv');

% sort them into chronological order.
rem_idx = zeros(1,length(d));
for iF = length(d):-1:1
    parts = strsplit(d(iF).name, {'_', '.'});
    
    inter_ver(iF) = str2double(parts{1});
    file_list{iF} = d(iF).name;
    
    f_num(iF) = str2double( regexp(d(iF).name, '\d+','match', 'once' )); % used to sort the files.
end

[~, idx] = sort(f_num);              % sort the file names based on actual numbers not first string value.
file_list  = file_list(idx)'

%% read all the files and save the

% set up the import options as per the ezTrack output format.
opts = detectImportOptions(file_list{1},'VariableNamesLine',1);
opts.DataLines = [2 inf]; % correct for one header row.
% find the variables of interest
x_idx = find(contains(opts.SelectedVariableNames, 'X'));
y_idx = find(contains(opts.SelectedVariableNames, 'Y'));
dist_idx = find(contains(opts.SelectedVariableNames, 'Distance_cm'));

if ~isempty(dist_idx)
    opts.SelectedVariableNames = ([x_idx, y_idx, dist_idx]); % should correspond to 'frames', 'X', 'Y', 'Distance_px'
else
    opts.SelectedVariableNames = ([x_idx, y_idx]); % should correspond to 'frames', 'X', 'Y', 'Distance_px'
end
opts.Delimiter = {','}; % required or else it assumes it is a space delimiter.
opts = setvartype(opts, 'double');

% sep up some varibales to be filled
x_y = [];
dist = [];
file_n = []; % video number
for ii = 1:length(file_list) % import the data
    t_table = table2array(readtable(file_list{ii},opts));
    x_y = [x_y ; t_table(:,1:2)]; % collect x y position
    if ~isempty(dist_idx)
        dist = [dist ; t_table(:,3)];
    end
    parts = strsplit(file_list{ii}, '_');
    file_n = [file_n; repmat(str2double(parts{1}), size(t_table,1),1)];
end

% same thing but for the ROI which is a string.

opts = detectImportOptions(file_list{1},'VariableNamesLine',1);
opts.DataLines = [2 inf]; % correct for one header row.
roi_idx = find(contains(opts.SelectedVariableNames, 'ROI_label'));
if ~isempty(roi_idx)
    opts.SelectedVariableNames = (roi_idx); % should correspond to 'ROI_label'
    opts.Delimiter = {','}; % required or else it assumes it is a space delimiter.
    opts = setvartype(opts, 'string');
    % sep up some varibales to be filled
    ROI = {};
    
    for ii = 1:length(file_list)
        t_table = table2array (readtable(file_list{ii},opts));
        ROI = [ROI ; t_table(:,1)];
        
    end
else
    ROI = [];
end
% check that the files are in the right order
if sum(unique(diff(file_n)) > 1) >0
    warning('Video order may not be correct')
end

%% load the timestamp file and fill in the tvec

TS = readtable('timeStamps.csv');

tvec = table2array(TS(:,2));

% correct for offsets if needed
% if tvec(1) ~= 0
%     tvec = tvec+abs(tvec(1));
% end
tvec = tvec./1000; % convert to seconds
%% convert to pos TSD

pos = tsd(tvec, x_y', {'x', 'y'});
pos.units = 'px';

% if distance is in cm, then convert distance to speed and smooth. 
if ~isempty(dist_idx)
    velo = dist./mode(diff(tvec)); 

    velo_smooth = smoothdata(velo, 'gaussian', round(1/mode(diff(tvec))));
    
    spd = tsd(tvec, [velo, velo_smooth]', {'Speed', 'Speed_smooth_win_1s'});
else
    spd = []; 
end
%% check plot

if plot_flag
    c_ord = winter(length(pos.tvec));
    subplot(2,2,[1 3])
    scatter(pos.data(1,:), pos.data(2,:), repmat(5, size(pos.tvec)), c_ord)
    colormap(c_ord)
    c_ax = colorbar(gca, 'Location', 'northoutside', 'Orientation', 'horizontal');
    c_ax.Label.String = 'time';
    c_ax.Ticks = [0 1];
    c_ax.TickLabels = {'Start' 'End'};
    subplot(2,2,2)
    plot(pos.tvec, pos.data)
    ylabel('position')
    legend('x', 'y')
    if ~isempty(dist_idx)
        subplot(2,2,4)
        plot(pos.tvec,spd.data(2,:))
        ylabel('Speed (cm/s)')
        xlabel('time (s)')
    end
    SetFigure([], gcf)
    maximize
    saveas(gcf, 'ezTrack_output.png')
    saveas(gcf, 'ezTrack_output.fig')
    
end
