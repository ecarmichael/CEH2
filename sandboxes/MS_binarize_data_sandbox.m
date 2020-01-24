function ms_data_in = MS_binarize_data_sandbox(cfg_in, ms_data_in)
%% MS_binarize_data: convert the ms_data into a binary (1 when above threshold, 0 if below). 
%
%
%
%    Inputs: 
%     - cfg_in [struct] contains user configurations (see below for
%     defaults)
%     - ms_data_in: [struct] this is the 'ms' output from the Miniscope
%     loader pipeline.  Should contain raw and filtered traces.  
%
%
%
%    Outputs: 
%     - data_out [nCell x time]: binary version of the input data.  
%
%   TODO
%     - for the rise find a way to detect multiple peaks that do not go
%     back down below the threshold. (might be possible in FindPeaks)
%
% EC 2020-01-21   initial version 
%
%
%% initialize

cfg_def = [];
cfg_def.method = 'zscore'; % can also be 'raw' to use a specific threshold value. 
cfg_def.rise = []; %Use 'cfg.rise = 1' to only use point where the data crosses the threshold. 
cfg_def.threshold = 2; % sd using zscore
cfg_def.operator = '>'; % which direction. 
cfg_def.data_type = 'RawTraces'; % can be RawTraces or FiltTraces, but can be others if they already exist in the ms structure. 
cfg_def.rise.mindistance = 300 ; % seems to work

cfg = ProcessConfig(cfg_def, cfg_in);

%% loop through cells and zscore the data. 

for iCell = size(ms_data_in.(cfg.data_type),2):-1:1
    
    switch cfg.method
        
        case 'zscore'
            data_temp(:,iCell) = zscore(ms_data_in.(cfg.data_type)(:,iCell));
            
        case 'raw'
            disp('Using RAW values') 
            data_temp(:,iCell) =ms_data_in.(cfg.data_type)(:,iCell);

        case 'rise'
            warning('WIP...');
            
%             [~, idx] = findpeaks(data_temp(:,iCell), 'minpeakdistance', cfg.rise.mindistance,'minpeakheight',  2); 
            
    end
    
end



%% convert to binary based on threshold 

switch cfg.operator
    
    case '>'
        data_temp(data_temp <= cfg.threshold) = 0;
        data_temp(data_temp > cfg.threshold) = 1;

        
    case '<'
        data_temp(data_temp >= cfg.threshold) = 0;
        data_temp(data_temp < cfg.threshold) = 1;
end

%% clean up

ms_data_in.BinaryTraces = data_temp;
% ms.hi
if isfield(ms_data_in, 'history')
        ms_data_in.history.fun_name{end+1} = {sprintf('MS_binarize_data on %s', date)};
        ms_data_in.history.cfg{end+1} = cfg; 
else
        ms_data_in.history.fun_name = {sprintf('MS_binarize_data on %s', date)};
        ms_data_in.history.cfg = {cfg}; 
end






