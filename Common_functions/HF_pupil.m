function [eye_tsd, blink_iv] = HF_pupil(dlc_dir, filt_flag, plot_flag)
%% HF_pupil: loads pupil position data from DLC and computes the diameter, direction (TBD), and blinks for the eye
%
%
%
%    Inputs: 
%    - dlc_dir: [path] location of the processed DLC videos and timestamps
%
%    - filt_flag: [logical] do you want to use the filtered predictions? 1
%    if yes, 0 if no, [no is default]
%
%
%
%    Outputs: 
%    - eye_tsd: [struct]  time-series data ('tsd') format. 
%                   eye_tsd.tvec: time
%                   eye_tsd.data(1,:): pupil diameter
%                   eye_tsd.data(2:3,:): gaze direction (x, y)
%    - blink_iv: blink times in the interval 'iv' format
%
%
%
% EC 2026-01-15   initial version 
%
%
%
%% initialize

if nargin < 1
    dlc_dir = cd; 
    filt_flag = 0;
    plot_flag = 0; 
elseif nargin < 2
    filt_flag = 0; 
    plot_flag = 0;
end


%% load the DLC data

eye = MS_DLC2TSD(dlc_dir, ['filtered'], [1 1], 0); 



%%  get the distance between opposite points


%% check the data

if plot_flag == 1


    figure(8919)
    clf
    hold on
    imagesc(eye.mean_frame)
    c_ord = MS_linspecer(ceil(size(eye.data,1)/2)); 
for ii = 1:10:1000
    cla
        imagesc(eye.mean_frame)

sc = scatter(eye.data(1:2:end,ii), eye.data(2:2:end,ii), 50, 'red', 'filled'); 
pause(.1)
clear sc
drawnow
end


end


