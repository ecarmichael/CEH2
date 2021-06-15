function [evts_out, IV_out] = MS_get_events(binary_in) 
%% MS_get_events: identify all binary blocks 
%
%
%
%    Inputs: 
%    - binary_in [1 x nSamples] binary.  
%
%    Outputs: 
%    - evts    [nEvents x 2]  start and stop indicies for each event
%    transient crossing the threshold. 
%
%    - IV [ struct]  this is the IV 'interval data' format from the mvdMlab
%    codebase.  allows use of the all powerful IV functions if using
%    vandermeerlab codebase. 
%
%
% EC 2021-06-15   initial version 
%   based off of MS_get_ca_events.  This function is simplified. 
%
%


%% find the peaks

[~,starts,~] = findpeaks(diff(binary_in), 'Threshold', .2); % get the start of the peak.  

[~,ends,~] = findpeaks(-diff(binary_in), 'Threshold', .2); % get the start of the peak.  


% correct for diff. 

if length(starts) < length(ends) 
    starts = [1; starts+1];
    start_flag = 1; 
else
    start_flag = 0; 
end

if length(starts) > length(ends) % correct for events at the end of the recording.  
    ends = [ends; length(binary_in)];
        endflag = 1; 
else
    endflag = 0; 
end

% put them together. 
    evts_out = [starts, ends]; 


IV_out.type = 'iv';
IV_out.tstart = starts;
IV_out.tend = ends; 
IV_out.usr = [];
% log this function. 
IV_out.cfg.history.mfun{1} = mfilename;
IV_out.cfg.history.cfg{1} = [];
%% plot for debugging

figure
hold on
plot(binary_in)
plot(evts_out(:,1), ones(1,length(evts_out))+.1, 'xr')
plot(evts_out(:,2), ones(1,length(evts_out))+.1, 'ob')
