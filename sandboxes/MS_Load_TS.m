function TS_out = MS_Load_TS(cfg_in)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   MS_Load_TS: imports timestamp.dat files from MiniScope data.  This will
%   identify all the timestamp* files present in the current directory.  
%
%
%   To Do: use a case based system to adapt to how the timestamp files are
%   stored. At the moment it assumes they are in the same folder and labed
%   in order.  
%
%       option 1: run it as part of a function that goes into all the
%       folders and pulls them out. 
%
%       option 2: have the timestamp files be appended to the ms data
%       struct from the preprocessing pipeline. 
%
%       option 3: use a hard coded list of files to find (in order). 
%  
%
% Example:
%   [camNum,frameNum,sysClock,buffer1] = importfile('timestamp.dat',2, 2952);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2019/11/22 14:53:51

%% Initialize variables.

cfg_def = [];
cfg_def.fname = {};            % put the file names that you want in here: {'timestamp1.dat', 'timstamp2.dat'} EC: make this adaptive to other file paths
cfg_def.correct_time = 1;

cfg = ProcessConfig2(cfg_def, cfg_in);


delimiter = '\t';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%[^\n\r]';

%% loop over timestamp.dat files found in this dir


if isempty(cfg.fname)
    ts_files_here = dir('timestamp*');
    
    % deal with more than 9 TS files
    for ii = 1:length(ts_files_here)
        temp = regexp(ts_files_here(ii).name, '\d*','Match');
        temp_ids(ii) = str2num(temp{1});
    end
    
    [~, idx] = sort(temp_ids);              % sort the file names based on number values
    ts_files_here = ts_files_here(idx);
    
    ts_files_here = {ts_files_here.name}'; % pull out a list of file names. 
    
elseif isa(cfg.fname, 'cell')
    ts_files_here = cfg.fname;  % take the cell array of file names. The order you put in is the order of TS cells you get out. 
else
    error('cfg.fname is not a cell or cell array.  If you want all the *timestamp files in a dir leave cfg.fname empty or don''t include it')
end






for iFiles = 1:length(ts_files_here)


% Open the text file.
fileID = fopen(ts_files_here{iFiles},'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);


%% Allocate imported array to column variable names
camNum = dataArray{:, 1};
frameNum = dataArray{:, 2};
sysClock = dataArray{:, 3};
buffer1 = dataArray{:, 4};

%% Correct the timeing to the first data point

if cfg.correct_time
    sysClock = [sysClock(1); sysClock(2:end) + sysClock(1)];
    disp('Correcting timestamps to first value in file')
end

%% extract multiple camera streams and format in the 'ts' manner
% which timestamp file is this? of how many

TS_out{iFiles}.type = 'ts';
TS_out{iFiles}.unit = 'ms';
TS_out{iFiles}.filename = ts_files_here{iFiles};
TS_out{iFiles}.file_num = iFiles;
TS_out{iFiles}.file_total = length(ts_files_here);
if isempty(cfg.fname)
    TS{iFiles}.file_total_flag = 'used given input so total is only from that input array';
else
    TS{iFiles}.file_total_flag = 'Used all timestamp* files in dir to get total';
end
TS_out{iFiles}.cfg.history.mfun = 'MS_Load_TS';
TS_out{iFiles}.cfg.history.cfg = cfg;
TS_out{iFiles}.cfg.filename = fullfile(pwd,TS_out{iFiles}.filename);

this_cam = 0;
for ii = unique(camNum)'
    this_cam = this_cam+1;
    idx = camNum == ii;
    TS_out{iFiles}.camera_id{this_cam} = camNum(idx);
    TS_out{iFiles}.framenumber{this_cam} = frameNum(idx);
    TS_out{iFiles}.system_clock{this_cam} = sysClock(idx);
    TS_out{iFiles}.buffer{this_cam} = buffer1(idx);
    TS_out{iFiles}.cfg.Fs{this_cam} = 1/(median(diff(TS_out{iFiles}.system_clock{this_cam}(2:end)))*0.001);
end

end


