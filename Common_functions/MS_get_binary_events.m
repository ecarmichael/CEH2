function evts = MS_get_binary_events(data_in, endflag) 
%% MS_get_ca_events: identify all the Ca transients cross a threshold. 
%
%
%
%    Inputs: 
%    - data_in [1 x nSamples] binary signal.  
%
%    - endflag [logical]   (optional) : what to do with events that do not end before the
%    end of the recording.  can be 0 (default) or 1.  If 0 then set end of last event as last index. If 1 then exclude the
%    last event.
%
%
%    Outputs: 
%    - evts    [nEvents x 2]  start and stop indicies for each Ca
%    transient crossing the threshold. 
%
%
%
%
% EC 2020-12-18   initial version 
%
%
%
%% initialize

if nargin == 1
    endflag= 0; % 2 SD as the default.  
end



%% find the peaks

[~,starts,~] = findpeaks(diff(data_in), 'Threshold', .2); % get the start of the peak.  

[~,ends,~] = findpeaks(-diff(data_in), 'Threshold', .2); % get the start of the peak.  


% correct for diff. 
evts(:,1) = starts+1; 

if length(starts) > length(ends) % correct for events at the end of the recording.  
        evts(:,2) = [ends; length(data_in)];
else
    evts(:,2) = ends;
end

if endflag
    evts = evts(:,1:end-1);
end

%% plot for debugging


% figure
% hold on
% plot(data_in)
% plot(evts(:,1), ones(1,length(evts))+.1, 'xr')
% plot(evts(:,2), ones(1,length(evts))+.1, 'ob')


