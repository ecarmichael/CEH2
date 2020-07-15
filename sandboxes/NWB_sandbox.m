% NDW/B sandbox
nwb_dir = 'C:\Users\ecarm\Documents\GitHub\matnwb';
%%
cd(nwb_dir)
addpath(genpath(pwd));
generateCore(); % generate core namespace located in the repository.


%% load a csv file from the google sheet for each animal

% meta = xlsread('J:\Williams_Lab\II_classification\Sub_log - Sample.csv')
% 
% fid = fopen('J:\Williams_Lab\II_classification\Sub_log - Sample.csv', 'rt');  %the 't' is important!
% C = textscan('%f,%f','HeaderLines',8);
% fclose(fid);

opts = detectImportOptions('J:\Williams_Lab\II_classification\Sub_log - R153.csv');
meta_table = table2cell(readtable('J:\Williams_Lab\II_classification\Sub_log - R153.csv', opts));
meta = [];
val_name = meta_table(:,1); 
for iR = 1:length(meta_table)
    if strcmp(val_name{iR}, 'Surgery')
        meta.Surgery.AP = meta_table{iR,3};
        meta.Surgery.ML = meta_table{iR,4};
        meta.Surgery.DV = meta_table{iR,5};
    else
        meta.(strrep(val_name{iR}, ' ', '_')) = meta_table{iR,2}; 
    end
end
%%
nwb = NwbFile( ...
    'session_description', meta.Task,...
    'identifier', meta.Subject_ID, ...
    'session_start_time', datetime([datestr(meta.Date, 'yyyy-mm-dd') ' ' datestr(meta.Start_time, 'HH:MM:ss')],'InputFormat', 'yyyy-MM-dd HH:mm:ss'), ...
    'general_experimenter', meta.Experimenter, ... % optional
    'general_session_id', num2str(meta.Session), ... % optional
    'general_institution', meta.Institution, ... % optional
    'general_related_publications', ''); % optional

subject = types.core.Subject( ...
    'subject_id', meta.Subject_ID, ...
    'age', ['P' num2str(meta.Subject_Age) 'D'], ...
    'description', [meta.Subject_ID '_' strrep(meta.Task, ' ', '_')], ...
    'species', meta.Species, ...
    'sex', meta.Sex);

nwb.general_subject = subject; 

%% position data
position_data = [linspace(0,10,100); linspace(0,8,100)]';
spatial_series_ts = types.core.SpatialSeries( ...
    'data', position_data, ...
    'reference_frame', '(0,0) is bottom left corner', ...
    'timestamps', linspace(0, 100)/200); 

