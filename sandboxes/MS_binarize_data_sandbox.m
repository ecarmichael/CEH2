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
cfg_def.threshold = 2; % sd using zscore
cfg_def.operator = '>'; % which direction.
cfg_def.data_type = 'RawTraces'; % can be RawTraces or FiltTraces, but can be others if they already exist in the ms structure.

% for the rise detection
cfg_def.rise.check_cell = 1;
cfg_def.rise.check = 1; % toggle check plot;
cfg_def.rise.smooth_type = 'sgolay'; % can be 'gaussian' (uses cfg.rise.gaussian configs) OR 'sgolay' ["Savitzky-Golay (polynomial) smoothing filter"] using cfg.rise.sgolay configs
cfg_def.rise.gaussian.smooth_factor = 10;
cfg_def.rise.mindistance = 200 ; % seems to work
cfg_def.rise.sgolay.order = 8; % best between 1-11;
cfg_def.rise.sgolay.framenum = 35; % best between 9-45

rise = ProcessConfig(cfg_def.rise, cfg_in.rise); % get the rise components. if not done here they can be overwritten by an incomplete .rise field in the cfg_in;

cfg = ProcessConfig2(cfg_def, cfg_in);
cfg.rise = rise; % put the rise configs back in.
clear rise
%% loop through cells and zscore the data.
ts_out= ts(size(ms_data_in.(cfg.data_type),2)); % keep an empty ts struct to append each cell

for iCell = size(ms_data_in.(cfg.data_type),2):-1:1
    
    switch cfg.method
        
        case 'zscore'
            data_temp(:,iCell) = zscore(ms_data_in.(cfg.data_type)(:,iCell));

           % convert to binary based on threshold
            switch cfg.operator
                case '>'
                    data_temp(data_temp <= cfg.threshold) = 0;
                    data_temp(data_temp > cfg.threshold) = 1;
                case '<'
                    data_temp(data_temp >= cfg.threshold) = 0;
                    data_temp(data_temp < cfg.threshold) = 1;
            end
            
        case 'raw'
            disp('Using RAW values')
            data_temp(:,iCell) =ms_data_in.(cfg.data_type)(:,iCell);
           
            % convert to binary based on threshold
            switch cfg.operator
                case '>'
                    data_temp(data_temp <= cfg.threshold) = 0;
                    data_temp(data_temp > cfg.threshold) = 1;
                case '<'
                    data_temp(data_temp >= cfg.threshold) = 0;
                    data_temp(data_temp < cfg.threshold) = 1;
            end
            
        case 'rise'
            %             warning('WIP...');
            switch cfg.rise.smooth_type
                case 'none'
                    this_data_in = ms_data_in.(cfg.data_type)(:,iCell);
                    
                case 'gaussian'
                    this_data_in = smoothdata(ms_data_in.(cfg.data_type)(:,iCell),cfg.rise.smooth_type,cfg.rise.gaussian.smooth_factor);
                    
                case 'sgolay'
                    this_data_in =  sgolayfilt(ms_data_in.(cfg.data_type)(:,iCell),cfg.rise.sgolay.order, cfg.rise.sgolay.framenum);
                    
            end
            this_data_z_diff = zscore(diff(this_data_in)); % convert to diff to get peaks and zscore for thresholding.
            
            % get the binary version of the filtered data
            switch cfg.operator
                
                case '>'
                    keep_idx = this_data_z_diff > cfg.threshold;
                case '>='
                    keep_idx = this_data_z_diff >= cfg.threshold;
                case '<'
                    keep_idx = this_data_z_diff < cfg.threshold;
                case '<='
                    keep_idx = this_data_z_diff >= cfg.threshold;
                otherwise
                    error('Unknown cfg.operator. Should be either ''>'', ''>='', ''<'', or ''<=''')
            end
            % assign a binary value
            this_data_binary = this_data_z_diff;
            this_data_binary(keep_idx) = 1;
            this_data_binary(~keep_idx) = 0;
            
            % find the peaks in the binarized trace. (essentiall just
            % finiding the start of each block of Ca activity above the
            % threshold.
            [pks, idx] = findpeaks(zscore(this_data_binary), 'minpeakdistance', cfg.rise.mindistance,'minpeakheight',  cfg.threshold);
            
            
            if cfg.rise.check && iCell == cfg.rise.check_cell% check that the find peaks is more or less working.
                figure(111)
                c_ord = linspecer(3);
                subplot(2,1,1)
                hold on
                findpeaks(zscore(this_data_binary), 'minpeakdistance', cfg.rise.mindistance,'minpeakheight',  cfg.threshold)
                %                 h(1) = gca;
                h(2) =plot( this_data_in, 'color', c_ord(1,:)); % add in the time and correct for diffs when needed.
                h(3) = plot( this_data_z_diff, 'color', c_ord(2,:), 'linewidth', 1.5);
                h(4) = plot( this_data_binary, 'color', c_ord(3,:));
                h(5) = plot(idx,ones(1,length(idx))*cfg.threshold, '+', 'MarkerSize', 10);
                
                legend(h([3,2,4,5]),{['Filter: ' cfg.rise.smooth_type],['Binary of zscore + diff +' cfg.rise.smooth_type], ['zscore + diff +' cfg.rise.smooth_type], ['peak detect -diff+'  cfg.rise.smooth_type]});
                ylim([min(this_data_z_diff), max(this_data_z_diff)]);
                
                subplot(2,1,2)
                hold on
                findpeaks(zscore(this_data_binary), 'minpeakdistance', cfg.rise.mindistance,'minpeakheight',  cfg.threshold)
                %                 h(1) = gca;
                h(2) =plot( this_data_in, 'color', c_ord(1,:)); % add in the time and correct for diffs when needed.
                h(3) = plot( this_data_z_diff, 'color', c_ord(2,:));
                h(4) = plot( this_data_binary, 'color', c_ord(3,:));
                h(5) = plot(idx,ones(1,length(idx))*cfg.threshold, '+', 'MarkerSize', 10);
                
                legend(h([3,2,4,5]),{['Filter: ' cfg.rise.smooth_type],['Binary of zscore + diff +' cfg.rise.smooth_type], ['zscore + diff +' cfg.rise.smooth_type], ['peak detect -diff+'  cfg.rise.smooth_type]});
                xlim([(idx(end)-250), (idx(end) + 250)])
                ylim([-pks(end), pks(end)]);
                
                title(['Cell #' num2str(cfg.rise.check_cell) '  evt# ' num2str(length(idx))])
                
                Resize_figure(gcf, 1.5)
                %                 pause(2)
                %                 close
            end % end check.
            
            % output the transitions in the binarized signal as a binary
            data_temp(:, iCell) = [this_data_binary; 0];
            
            data_temp(idx, iCell) = 1;
            data_temp(~idx, iCell) =0;
            
            % Use the TS format to track the pseudo spikes; Helpful for later
            % make sure you put the cells back in the right order.
            this_ts = ts({ms_data_in.time(idx)});
            this_ts.label = {num2str(iCell)};
            
            ts_out.t{iCell} = this_ts.t{1};
            ts_out.label{iCell} = this_ts.label{1};
    end
    
end





%% clean up

ms_data_in.BinaryTraces = data_temp;

if strcmp(cfg.method, 'rise')
    ms_data_in.spike_ts = ts_out;
end
% ms.hi
if isfield(ms_data_in, 'history')
    ms_data_in.history.fun_name{end+1} = {sprintf('MS_binarize_data on %s', date)};
    ms_data_in.history.cfg{end+1} = cfg;
else
    ms_data_in.history.fun_name = {sprintf('MS_binarize_data on %s', date)};
    ms_data_in.history.cfg = {cfg};
end






