function Lap_screener(behav, ms, z_thresh)
%% Lap_screener: makes a simple plot of the Ca traces, position, speed, SPF, and track activity
%
%
%
%    Inputs: 
%    - behav   [struct] 
%
%    - ms       [struct]
%
%    Outputs: 
%    -
%
%
%
%
% EC 2020-08-04   initial version 
%
%
%
%% initialize

if nargin <2
    error('Requires behav and ms inputs')
elseif nargin ==2
    fprintf('<strong>%s:</strong> no z-score specified, using default 2SD\n', mfilename)
elseif nargin == 3    
    fprintf('<strong>%s:</strong> Using user specified zscore for binarizing %d \n', mfilename, z_thresh)
end


%% make a simple plot. 

frame_time = 0.001*(mode(diff((ms.time))));  % convert time frame to seconds.

left_idx = MS_get_direction(behav.position(:,1), -.1); % use -threshold for leftbound and + for right.
right_idx = MS_get_direction(behav.position(:,1), .1); % use -threshold for leftbound and + for right.


movement_thresh = 2.5; % in cm/s
movement_idx = behav.speed >movement_thresh; % get times when the animal was moving.
left_idx = left_idx & movement_idx; % only keep the indices when they are moving to the left.
right_idx = right_idx & movement_idx;

[L_laps, L_lap_start_idx, L_lap_end_idx] = MS_get_laps(left_idx, floor(1.5*(1/frame_time)),floor(10*(1/frame_time)));
[R_laps, R_lap_start_idx, R_lap_end_idx] = MS_get_laps(right_idx, floor(1.5*(1/frame_time)),floor(10*(1/frame_time)));