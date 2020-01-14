function [iv_out, durations] = MS_extract_NLX_blocks_sandbox(cfg_in, evt)
%% MS_extract_NLX_blocks: extracts peaks in the Neuralynx event files.  
%   Designed to separate discontinuous recordings. Mainly used for aligning
%   csc recordings to miniscope data using a TTL on the miniscope DAQ. 
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
%
%
%
%    Outputs: 
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
cfg_def.check = 0; % use 1 to toggle plotting the output. 
cfg = ProcessConfig2(cfg_def, cfg_in); 

%% identify peaks in  diff(evt.t{5}) marking transitions in the camera TTLs
peak_threshold =  (mean(diff(evt.t{cfg.t_chan}) +cfg.peak_threshold*std(diff(evt.t{cfg.t_chan}))));
[peaks, idx] = findpeaks(diff(evt.t{cfg.t_chan}), 'minpeakheight',peak_threshold, 'minpeakdistance', cfg.min_dist);

% add on to each end.  
if idx(1) > 5 % arbitrary at this point, but there should be a peak right at the beginging unless the nlx and the miniscope TTL are the same. 
idx = [1 idx];
peaks = [diff(evt.t{cfg.t_chan}(1:2)) peaks];
disp('No jump in time at the start.  Making the start of the TTLs the start of the first recording block')
end

fprintf('\nDetected Sampling Frequency mode = %0.3fHz , mean = %0.3fHz\n', mode(diff(evt.t{cfg.t_chan})),mode(diff(evt.t{cfg.t_chan})))
fprintf(['Detected %.0f trigger transitions treating this as %.0f distinct recordings\n'], length(idx), length(idx))

    x_val = evt.t{1}:mode(diff(evt.t{cfg.t_chan})):evt.t{2};


% convert to iv
tstart = evt.t{cfg.t_chan}(idx+1);
tend = [evt.t{cfg.t_chan}(idx(2:end))  evt.t{cfg.t_chan}(end)];

iv_out = iv(tstart, tend);


    durations = tend - tstart;
% durations = (tend - tstart).*mode(diff(evt.t{cfg.t_chan}));


% plot the diff and the detected peaks as a check.
if cfg.check == 1
    figure(1234)
    hold on
    subplot(2,1,1)
    plot(diff(evt.t{cfg.t_chan}), 'k')
    hline(peak_threshold, '--r')
    % plot(Rec_idx, 100, '*k')
    for iRec = 1:length(tstart)
        % text(Rec_idx(iRec),Rec_peak(iRec),num2str(iRec))
        text(tstart(iRec),peak_threshold,num2str(iRec))
        
    end
    
    c_ord = linspecer(length(idx)); % nice colours. 
    subplot(2,1,2) % was going to make a coole figure with all the events as lines on different y points across the entrie session. 
    evt_vals = evt.t{cfg.t_chan};
    plot(x_val, NaN(size(x_val)))
    hold on
    for iRec = 1:length(tstart)
            line([tstart(iRec), tend(iRec)], [iRec iRec],'color', c_ord(iRec,:), 'linewidth', 3)
        text(median([tstart(iRec), tend(iRec)]), iRec+0.5, sprintf('%2.1fs', durations(iRec)))
    end
    ylim([0 length(idx)+1]);
end

