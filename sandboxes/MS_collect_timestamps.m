function [ts_out, ts_fnames] = MS_collect_timestamps(dir_name, output_type);
%% MS_collect_timestamps: goes through all folders in a dir to find timestamp.dat files and collects them.
%  Output 'ts_out' can be in cells (ts_out{1:N}) or as a structure with the
%  file names. Cells is the default.
%
%
%   inputs:
%       -dir_name: [string] path to folders containing timestamp.dat files.
%       Typically this is the output folder from Miniscope recordings.
%       - output_type: [string]  can be either 'cells' or 'struct'
%
%   Outputs:
%       - ts_out: [1xnCell] or [struct] contains all the timestamp variables.
%       
%       -ts_fname: [1xnCell] contains the names of the folders where the
%       corresponding timestamp.dat file was found. 
%
%
%   EC 2020-01-13 : initial version
%
%

%% Defaults

if nargin == 0
    error('Must specify the directory that you wish to search')
elseif nargin <2
    output_type = 'cells';
    disp('No output_type specified. Using ''cells''');
end

    
ts_out = [];

%% get all the timestamp files in the dir_name

    
TS_files = dir(fullfile(dir_name, '**', 'timestamp*.dat'));  
    

%% sort based on the date of the folder.
f_time = [];
for iF = length(TS_files):-1:1
    parts = strsplit(TS_files(iF).folder, filesep);
    
    Ts = strsplit(parts{end}, '_'); 

    f_time{iF} = datestr([Ts{1}(2:end) ':' Ts{2}(2:end) ':' Ts{3}(2:end)], 'hh:MM:ss');

    fprintf('%s  - %s\n', parts{end}, f_time{iF} )
end

% [~, idx] = sort([TS_files.datenum]);              % sort the file names based on number values
[~, idx] = sort([f_time]);              % sort the file names based on 
TS_files  = TS_files(idx);

%% extract the values.
this_dir = pwd;

for iF  = 1:length(TS_files)
    fldr = strsplit(TS_files(iF).folder, filesep);
    cd([dir_name filesep fldr{end}])
    fprintf('\nExtracting %d/%d Timestamps.dat from %s...\n', iF, length(TS_files), fldr{end})
    
    cfg_ts = [];
    cfg_ts.correct_time = 1;
    
    switch output_type
        
        case 'cells' 
            temp = MS_Load_TS(cfg_ts);
            ts_out{iF} = temp{1}; 
            
        case 'Cells'
            temp = MS_Load_TS(cfg_ts);
            ts_out{iF} = temp{1}; 
            
        case 'struct' 
            temp = MS_Load_TS(cfg_ts);
            ts_out.(fldr{end}) = temp{1}; 

        case  'Struct'
            temp = MS_Load_TS(cfg_ts);
            ts_out.(fldr{end}) = temp{1}; 
    end
    ts_fnames{iF} = fldr{end};
    temp = []; 
end

%% clean up

cd(this_dir)
fprintf('Done. Output is ''%s''\n', output_type)



