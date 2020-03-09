function [event_blocks, iv_out, durations] = MS_extract_NLX_blocks_sandbox(cfg_in, evt)
%% MS_extract_NLX_blocks: extracts peaks in the Neuralynx event files.  
%   Designed to separate discontinuous recordings. Mainly used for aligning
%   csc recordings to miniscope data using a TTL on the miniscope DAQ. This
%   will then apply a gitter to isolate the event blocks with no major
%   jumps in sampling which can occur when the TTL is not working properly.
%   The gitter correction can be toggled with cfg.gitter = 0 or 1. 
%
%
%
%    Inputs: 
%     - cfg_in: [struct] contains user paramters
%           -[defaults]
%               - cfg.t_chan = 3; NLX uses 1 for starting times, 2 for
%               stopping, and 3 - n for TTL events. 
%               - cfg.peak_threshold = 0.05 % number of standard deviations from the
%               mean difference in the samples for this t_chan.
%               - cfg.min_dist = 10: minimum number of samples between
%               detected peaks.  
%               - cfg.check = 0: toggles the plot to check the detection. 
%               - cfg.gitter = 1:  applies a gitter to find the indices
%               without jumps in the sampling which occur with
%               discontinuous recording. 
%
%
%
%    Outputs: 
%     - event_blocks: NLX events restricted to each detected segment. 
%
%     - iv_out: [struct] contains start and end times for each identified
%     block. 
%
%     - durations: [1 x N array] contains the duration of the event. 
%
%
%
%
% EC 2020-01-14   initial version 
%   Requies vandermeerlab codebase : https://github.com/vandermeerlab
%
%% set defaults
cfg_def = [];
cfg_def.t_chan = 3; % which channel in the evt.t to use.
cfg_def.peak_threshold = 0.05; % number of Std above the mean of the diff(evts)
cfg_def.min_dist = 10; % distance between peaks
cfg_def.check = 1; % use 1 to toggle plotting the output. 
% gitter configs
cfg_def.gitter = 1; % toggle use of gitter correction
cfg_def.gitter_threshol = 1.05; 
cfg_def.low_val = 10; % how many samples back you want to start the gitter. Making this too big could overlap with other segments
cfg_def.high_val = 10; % in samples. If this is too big it could overlap with subsequent events.
cfg_def.bad_block =[]; % known recording blocks to skip.  
cfg = ProcessConfig2(cfg_def, cfg_in); 

if ~isempty(cfg.bad_block)
    for ii = cfg.bad_block
        fprintf('<strong>%s</strong>: block: <strong>%d</strong> flagged for skipping.  Only using first pass detection, no gitter correction. Removed prior to resizing.\n', mfilename,ii)
    end
end
%% identify peaks in  diff(evt.t{5}) marking transitions in the camera TTLs
peak_threshold =  (mean(diff(evt.t{cfg.t_chan}) +cfg.peak_threshold*std(diff(evt.t{cfg.t_chan}))));
[peaks, idx] = findpeaks(abs(diff(evt.t{cfg.t_chan})), 'minpeakheight',peak_threshold, 'minpeakdistance', cfg.min_dist);

% add on to each end.  
if idx(1) > 5 % arbitrary at this point, but there should be a peak right at the beginging unless the nlx and the miniscope TTL are the same. 
idx = [1 idx];
peaks = [diff(evt.t{cfg.t_chan}(1:2)) peaks];
disp('No jump in time at the start.  Making the start of the TTLs the start of the first recording block')

tstart = evt.t{cfg.t_chan}(idx+1);
tstart(1) = evt.t{cfg.t_chan}(idx(1)); % correct for first timestamp. 

else
    tstart = evt.t{cfg.t_chan}(idx+1);
end



tend = [evt.t{cfg.t_chan}(idx(2:end))  evt.t{cfg.t_chan}(end)];
% convert to IV format. 
iv_out = iv(tstart, tend);
% get duration of each block. 
durations = tend - tstart;
% durations = (tend - tstart).*mode(diff(evt.t{cfg.t_chan}));

% print something about the detected event breaks. 
fprintf('\nDetected Sampling Frequency mode = %0.3fHz , mean = %0.3fHz\n', mode(diff(evt.t{cfg.t_chan})),mode(diff(evt.t{cfg.t_chan})))
fprintf(['Detected %.0f trigger transitions treating this as %.0f distinct recordings\n'], length(idx), length(idx))



% plot the diff and the detected peaks as a check.
if cfg.check == 1
    figure(1234)
    hold on
    subplot(5,1,1)
    plot(diff(evt.t{cfg.t_chan}), 'k')
    hline(peak_threshold, '--r')
    % plot(Rec_idx, 100, '*k')
    for iRec = 1:length(tstart)
        text(idx(iRec),peaks(iRec),num2str(iRec))
        %         text(tstart(iRec),peak_threshold,num2str(iRec))
    end
    title('Peak Detection')
    
    c_ord = linspecer(length(idx)); % nice colours.
    subplot(5,1,2:3) % was going to make a coole figure with all the events as lines on different y points across the entrie session.
    
    x_val = evt.t{1}:mode(diff(evt.t{cfg.t_chan})):evt.t{2};
    plot(x_val, NaN(size(x_val)))
    hold on
    for iRec = 1:length(tstart)
        line([tstart(iRec), tend(iRec)], [iRec iRec],'color', c_ord(iRec,:), 'linewidth', 3)
        text(tstart(iRec), iRec+0.5, sprintf('#%0.0f %2.1fs',iRec, durations(iRec)))
    end
    ylim([0 length(idx)+1]);
    title('Pre check')

end

%% apply a gitter to resolve jumps in timing. This entire section could be replaced with proper use 
%  of the max(diff(evt.t)) and using indexing, but this is safer since it can deal with multiple 
%  peaks. 

if cfg.gitter ==1
    disp('MS_extract_NLX_blocks_sandbox: Running gitter correction...')
    for iRec = 1:length(tstart)
    
        if ~ismember(iRec, cfg.bad_block)
            
%         if iRec == 1 && idx(iRec) ==1
%             temp_evt = restrict(evt, evt.t{cfg.t_chan}(idx(iRec)), evt.t{cfg.t_chan}(idx(iRec+1)));
%             nSamples_high = length(temp_evt.t{cfg.t_chan});
%  
%         elseif iRec == length(idx)
%             temp_evt = restrict(evt, evt.t{cfg.t_chan}(idx(iRec)), evt.t{cfg.t_chan}(end));
%             nSamples_low = length(temp_evt.t{cfg.t_chan});
% 
%         else
%             temp_evt = restrict(evt, evt.t{cfg.t_chan}(idx(iRec)), evt.t{cfg.t_chan}(idx(iRec+1)));
%             nSamples = length(temp_evt.t{cfg.t_chan});
%         end
%         
%         if nSamples < (cfg.low_val + cfg.high_val) 
%             low_val = 0; 
%             high_val = 0;
%         else
            low_val = cfg.low_val;
            high_val = cfg.high_val;
%         end
    idx_low = NaN; % initialize with something.
    idx_high = NaN;
    
    if iRec == 1 && idx(iRec) ==1
        low_val = 0; 
        while ~isempty(idx_low)
            low_val = low_val - 1;
            temp_evt = restrict(evt, evt.t{cfg.t_chan}(idx(iRec)-low_val), evt.t{cfg.t_chan}(idx(iRec+1)));
            idx_low =find(diff(temp_evt.t{cfg.t_chan}(1:20)) > mode(diff(temp_evt.t{cfg.t_chan}))*cfg.gitter_threshol);
        end
        
        
        fprintf('\nMS_extract_NLX_blocks_sandbox: First block starts at the begining of the recoring (idx = %.0f).\n', low_val)
   
         while ~isempty(idx_high)
            high_val = high_val -1;
            temp_evt = restrict(evt, evt.t{cfg.t_chan}(idx(iRec)), evt.t{cfg.t_chan}(idx(iRec+1)+high_val));
            idx_high =find(diff(temp_evt.t{cfg.t_chan}(end-20:end)) > mode(diff(temp_evt.t{cfg.t_chan}))*cfg.gitter_threshol);
            
%             figure(11)
%             plot(diff(temp_evt.t{cfg.t_chan}))
%             ylim([0 mode(diff(temp_evt.t{cfg.t_chan}))+1])
%             pause(.5)
%             close
        end
        
        
    elseif iRec == 1 && (idx(iRec) ~=1 && idx(iRec) < cfg.low_val) % use this block if the first index is within cfg.low_val of the start of the events.   
        low_val = (idx(iRec)-2);  % start as far back as the first spot will allow up to the cfg.low_val

        
         while ~isempty(idx_low)
            low_val = low_val - 1;
            temp_evt = restrict(evt, evt.t{cfg.t_chan}(idx(iRec)-low_val), evt.t{cfg.t_chan}(idx(iRec+1)));
            idx_low =find(diff(temp_evt.t{cfg.t_chan}(1:20)) > mode(diff(temp_evt.t{cfg.t_chan}))*cfg.gitter_threshol);
        end
        
         while ~isempty(idx_high)
            high_val = high_val -1;
            temp_evt = restrict(evt, evt.t{cfg.t_chan}(idx(iRec)-low_val), evt.t{cfg.t_chan}(idx(iRec+1)+high_val));
            idx_high =find(diff(temp_evt.t{cfg.t_chan}(end-20:end)) > mode(diff(temp_evt.t{cfg.t_chan}))*cfg.gitter_threshol);
        end
          
        
        
    elseif iRec == length(idx)
        high_val =1; % only move backwards from the end. this is 1 not zero since we start by subtracting 1 in each loop

        while ~isempty(idx_low)
            low_val = low_val - 1;
            temp_evt = restrict(evt, evt.t{cfg.t_chan}(idx(iRec)-low_val), evt.t{cfg.t_chan}(end));
            idx_low =find(diff(temp_evt.t{cfg.t_chan}(1:20)) > mode(diff(temp_evt.t{cfg.t_chan}))*cfg.gitter_threshol);
        end
        
         while ~isempty(idx_high)
            high_val = high_val -1;
            temp_evt = restrict(evt, evt.t{cfg.t_chan}(idx(iRec)-low_val), evt.t{cfg.t_chan}(end+high_val));
            idx_high =find(diff(temp_evt.t{cfg.t_chan}(end-20:end)) > mode(diff(temp_evt.t{cfg.t_chan}))*cfg.gitter_threshol);
        end
          
    else
        % loop until you find the index offset that gives no jump in time. for
        % start
         while ~isempty(idx_low)
            low_val = low_val - 1;
            temp_evt = restrict(evt, evt.t{cfg.t_chan}(idx(iRec)-low_val), evt.t{cfg.t_chan}(idx(iRec+1)));
            idx_low =find(diff(temp_evt.t{cfg.t_chan}(1:20)) > mode(diff(temp_evt.t{cfg.t_chan}))*cfg.gitter_threshol);

            % check  plot
%             figure(11)
%             plot(diff(temp_evt.t{cfg.t_chan}))
%             ylim([0 mode(diff(temp_evt.t{cfg.t_chan}))+1])
%             pause(.5)
%             close
         end
        
        % loop until you find the index offset that gives no jump in time. for
        % end
        while ~isempty(idx_high)
            high_val = high_val -1;
            temp_evt = restrict(evt, evt.t{cfg.t_chan}(idx(iRec)-low_val), evt.t{cfg.t_chan}(idx(iRec+1)+high_val));
            idx_high =find(diff(temp_evt.t{cfg.t_chan}(end-20:end)) > mode(diff(temp_evt.t{cfg.t_chan}))*cfg.gitter_threshol);
%             figure(11)
%             plot(diff(temp_evt.t{cfg.t_chan}))
%             ylim([0 mode(diff(temp_evt.t{cfg.t_chan}))+1])
%             pause(.5)
%             close
        end
        
%         disp(['Corrected indexting ' num2str(low_val) ' - ' num2str(high_val)])
    end % end iRec first/last/others
    
        else
            low_val = 1; 
            high_val = 1; 
        end
    times_to_use(iRec,:) = [low_val, high_val]; % keep the good values for restrictions in the next cell.
    
    if iRec == length(idx)
        temp_evt = restrict(evt, evt.t{cfg.t_chan}(idx(iRec)-low_val), evt.t{cfg.t_chan}(length(evt.t{cfg.t_chan}))+high_val);
    else
        temp_evt = restrict(evt, evt.t{cfg.t_chan}(idx(iRec)-low_val), evt.t{cfg.t_chan}(idx(iRec+1)+high_val));
    end
    biggest_jumps(iRec) = max(diff(temp_evt.t{cfg.t_chan}));

    fprintf('Evt: %.0f Best offsets: start = %.0f  end = %.0f largest jump = %.3f sec\n', iRec, times_to_use(iRec,1), times_to_use(iRec, 2), biggest_jumps(iRec))
    
    event_blocks{iRec} = temp_evt; % hold the resticted events for output
        durations(iRec) = event_blocks{iRec}.t{cfg.t_chan}(end) - event_blocks{iRec}.t{cfg.t_chan}(1);

    end % end number of start/stop indices 'idx'
    
if cfg.check == 1
    figure(1234)
    subplot(5,1,4:5) % was going to make a coole figure with all the events as lines on different y points across the entrie session. 
    x_val = evt.t{1}:mode(diff(evt.t{cfg.t_chan})):evt.t{2};
    plot(x_val, NaN(size(x_val)))
    hold on
    for iRec = 1:length(idx)
            line([event_blocks{iRec}.t{cfg.t_chan}(1),event_blocks{iRec}.t{cfg.t_chan}(end)], [iRec iRec],'color', c_ord(iRec,:), 'linewidth', 3)
        text(event_blocks{iRec}.t{cfg.t_chan}(1), iRec+0.5, sprintf('#%0.0f %2.1fs',iRec, durations(iRec)))
    end
    ylim([0 length(idx)+1]);
    title('Post gitter correction')
    Resize_figure(gcf, 1.5); % dumb function that resizes to 1.5x
end

    
else % use tehe raw values from the peak detection *[NOT RECOMMENEDED]*
    for iRec = 1:length(idx)
        
    event_blocks{iRec} = restrict(evt,tstart(iRec), tend(iRec));
    warning('Not using gitter sample correction')
    end
end% end gitter section;


% %% clean up
% 
% if isfield(ms_data_in, 'history')
%         ms_data_in.history.fun_name{end+1} = {sprintf('MS_extract_NLX_blocks on %s', date)};
%         ms_data_in.history.cfg{end+1} = cfg; 
% else
%         ms_data_in.history.fun_name = {sprintf('MS_extract_NLX_blocksa on %s', date)};
%         ms_data_in.history.cfg = {cfg}; 
% end

