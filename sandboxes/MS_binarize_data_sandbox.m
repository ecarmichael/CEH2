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

% for the rise detection
cfg_def.rise.check_cell = 1;
cfg_def.rise.check = 1; % toggle check plot; 
cfg_def.rise.smooth_type = 'none'; % can be 'gaussian' (uses cfg.rise.gaussian configs) OR 'sgolay' ["Savitzky-Golay (polynomial) smoothing filter"] using cfg.rise.sgolay configs
cfg_def.rise.gaussian.smooth_factor = 10; 
cfg_def.rise.mindistance = 200 ; % seems to work
cfg_def.rise.sgolay.order = 8; % best between 1-11; 
cfg_def.rise.sgolay.framenum = 17; % best between 9-45 

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
            %             warning('WIP...');
            switch cfg.rise.smooth_type
                case 'none'
                   this_data = ms_data_in.(cfg.data_type)(:,iCell); 

                case 'gaussian'
                    this_data = smoothdata(ms_data_in.(cfg.data_type)(:,cfg.rise.check_cell),cfg.rise.smooth_type,cfg.rise.gaussian.smooth_factor);
                    
                case 'sgolay'
                    this_data =  sgolayfilt(ms_data_in.(cfg.data_type)(:,cfg.rise.check_cell),cfg.rise.sgolay.order, cfg.rise.sgolay.framenum);
                    
            end
            [~, idx] = findpeaks(diff(zscore(this_data)), 'minpeakdistance', cfg.rise.mindistance,'minpeakheight',  cfg.threshold);
            
                    
            if cfg.rise.check && iCell == cfg.rise.check_cell% check that the find peaks is more or less working.
                figure(111)
                subplot(2,1,1)
                findpeaks(zscore(diff(this_data(:,cfg.rise.check_cell))), 'minpeakdistance', cfg.rise.mindistance,'minpeakheight',  cfg.threshold)
                hold on
                plot(zscore(diff(this_data(:,cfg.rise.check_cell))))
%                 plot(zscore((this_data(:,cfg.rise.check_cell)*10)) - mean(zscore(this_data(:,cfg.rise.check_cell))))
                legend({'peak detect -diff+smooth','detected peak', 'smooth og'});
                
                subplot(2,1,2)
                findpeaks(diff(zscore(this_data(:,cfg.rise.check_cell))), 'minpeakdistance', cfg.rise.mindistance,'minpeakheight',  cfg.threshold)
                [pks, idx] = findpeaks(diff(zscore(this_data(:,cfg.rise.check_cell))), 'minpeakdistance', cfg.rise.mindistance,'minpeakheight',  cfg.threshold);
                hold on
                plot(zscore(this_data(:,cfg.rise.check_cell)))
%                 plot((zscore(this_data(:,cfg.rise.check_cell))*10) - mean(zscore(this_data(:,cfg.rise.check_cell))))
                legend({'peak detect -diff+smooth','detected peak', 'smooth og'});
                xlim([(idx(end)-250), (idx(end) + 250)])
                ylim([-pks(end), pks(end)]);
                plot(idx, cfg.threshold, '*')
                title(['Cell #' num2str(cfg.rise.check_cell)])
                %                 pause(2)
                %                 close
            end % end check.
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






