function evts = MS_get_binary_events(data_in) 
%% MS_get_ca_events: identify all the Ca transients cross a threshold. 
%
%
%
%    Inputs: 
%    - data_in [1 x nSamples] binary signal.  

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
    thresh = 2; % 2 SD as the default.  
end



%% find the peaks

[~,starts,~] = findpeaks(diff(data_in), 'Threshold', .2); % get the start of the peak.  

[~,ends,~] = findpeaks(-diff(data_in), 'Threshold', .2); % get the start of the peak.  


% correct for diff. 
evts(:,1) = starts+1; 
evts(:,2) = ends;


%% plot for debugging


% figure
% hold on
% plot(data_in)
% plot(evts(:,1), ones(1,length(evts))+.1, 'xr')
% plot(evts(:,2), ones(1,length(evts))+.1, 'ob')


