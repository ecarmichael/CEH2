function [hypno_labels, time_labels] = MS_get_hypno_label(cfg_in, TS_name)
%% MS_get_hypno_label: extracts the sleep state identifier for each timestamp file.  
%
%
%
%    Inputs: 
%     - cfg [struct] configuration parameters.  See below for defauls. 
%     - TS_name [1xncells]  cells containing the directory names where the TS
%     timestamps.dat was extracted as [TS, TS_name] = MS_collect_timestamps 
%
%
%
%    Outputs: 
%     - labels [1xn cells] cells with names which best match the 'REM',
%     'Awake', "Task', 'SW' labels.  these can also be user specified
%     depending on what you are looking for. 
%
%
%
%
% EC 2020-01-21   initial version 
%
%
%% initialize
cfg_def = [];
cfg_def.label_to_find = {'SW', 'REM', 'Other'};

cfg = ProcessConfig(cfg_def, cfg_in); 

%% find a label for each name
hypno_labels = cell(1, length(TS_name));
time_labels = cell(1,length(TS_name));
for iT = 1:length(TS_name)
    
    lab_out = [];
    for iL = 1:length(cfg.label_to_find)
        
        if contains(TS_name{iT}, cfg.label_to_find{iL}, 'IgnoreCase', true)
            if contains(TS_name{iT}, 'REM', 'IgnoreCase', true)
                rem_idx = strfind(lower(TS_name{iT}), lower('REM'));
%                 if length(TS_name{iT}) >= rem_idx+3  % why did I do this
%                 (EC)? 
%                     [~, status]= str2num(TS_name{iT}(rem_idx+3));
%                     if status ~= 1
%                         lab_out = 'Other';
%                     else
%                         lab_out = cfg.label_to_find{iL};
%                     end
%                 else
                    lab_out = cfg.label_to_find{iL};
%                 end
            else
                lab_out = cfg.label_to_find{iL};
            end
        end
    end
    if isempty(lab_out)
        lab_out = 'Other';
    end
    hypno_labels{iT} = lab_out;
    % get the time string
    time_str = regexp(TS_name{iT},'\d*','Match');
    time_labels{iT} = datestr(datestr([time_str{1} ':' time_str{2} ':' time_str{3}]), 'HH:MM:SS');
    
end



