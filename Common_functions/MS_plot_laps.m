function MS_plot_laps(behav, laps, spike_times)
%% MS_plot_laps: makes a plot with laps as rows on a plot.  Spike_times will be added as dots on the position plot. 
%
%
%
%    Inputs: 
%    - behav [ struct] the behaviour structure from the miniscope animal
%    camera.  Will only use time and position variables. 
%
%    - laps  [1 x nLaps]  array with the indicies of each lap number (as
%    0000000111111100000002222222000000003333333....)
%
%    [Optional]
%    - spike_times  [n x nSpikes]  array with 'spike' times.  
%    Outputs: 
%
%
%    - h  [plot handle]
%
% EC 2020-07-16   initial version 
%
%
%
%% initialize

if nargin <2
    error('requires lap array')
elseif nargin == 2
    spike_times = []; % just leave it empty
end

%% setup the plot
nLaps = unique(laps);
nLaps(nLaps==0) = []; % remove the zero
max_y =(ceil(max(behav.position(:,2))*0.1))*10; % get the next largest 10 in the y axes for spacing. 
offsets = 1+(nLaps*max_y);
hold on
for iL = 1:length(nLaps)
    this_lap_idx = find(laps==nLaps(iL)); 
    plot(behav.position(this_lap_idx,1), behav.position(this_lap_idx,2)+offsets(iL), 'color', [0.8 0.8 0.8]);
    tick_val(iL) = median(behav.position(this_lap_idx,2)+offsets(iL)); 
    
    if ~isempty(spike_times)
        spike_idx = this_lap_idx(ismember(this_lap_idx,spike_times)); 
            plot(behav.position(spike_idx,1), behav.position(spike_idx,2)+offsets(iL),'.', 'color', [0.9153    0.2816    0.2878]); % make spikes appear as red points. 
    end
    
    
    
end
% ylabel('lap #');

set(gca, 'ytick', tick_val([1, end]), 'yticklabels', nLaps([1, end]))
set(gca, 'xtick', []);

