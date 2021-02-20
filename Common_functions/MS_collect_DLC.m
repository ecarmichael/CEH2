function data_out = MS_collect_DLC(dir_in, model_in)
%% MS_collect_DLC: loads and collects all DLC files in a directory. Will skip over files without a number since DLC saves the interation number in the .csv
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
%    Outputs:
%    - DLC_array [nTime x nParts]   concatenation of all the DLC files
%
%
%
%
% EC 2021-02-17   initial version
%
%
%
%% initialize

if nargin == 0
    dir_in = cd; % just use current dir.
    model_in = [];
elseif nargin == 1
    model_in = [];
end


%% find all the files

%%%%% TO Do implement model_in catch %%%%%%%%


og_dir = dir_in;
cd(dir_in);

file_list = {};
d = dir('*.csv');
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
end
    
% remove empty cells
rem_idx(rem_idx == 0) = [];
file_list(rem_idx) = [];

newest_inter_ver = max(inter_ver);
% loop and find only the DLC versions that use the newest model (ie most
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
    tbl = readtable(file_list{iF});
    DLC_labels = tbl{1,2:end};
    fields = unique(DLC_labels); % get the parts
    
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
for iF  = 1:length(file_list)
    for iFields = 1:length(fields)
        data_out.(fields{iFields}) = [ data_out.(fields{iFields});  this_field{iF, iFields}];
    end
end

cd(og_dir); 
