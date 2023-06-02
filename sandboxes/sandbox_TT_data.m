%% sandbox TT sample data

% add the vandermeer lab codebase 
addpath(genpath('/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared'));
disp('vandermeerlab codebase added')
% add EC's codebase here: 
addpath(genpath('/home/williamslab/Documents/Github/CEH2'))

data_dir = '/home/williamslab/Desktop/For Travis - LFP test/ground'; % where is the data

%% quick PSD from a session

MS_Quick_psd; % run this to get a PSD from all the channels in this dir. 


%% load things
cd(data_dir)

cfg = [];
cfg.fc = {'CSC2.ncs' 'CSC8.ncs', 'CSC10.ncs', 'CSC14.ncs'}; % comment this to load all channels in this folder. 

cfg.desired_sampling_frequency = 2000; % helps with speed. 

csc = MS_LoadCSC(cfg); % load csc channels. 

art_idx = zscore(abs(csc.data)) > 3; % get index for large amp events, likely hitting a wall or something. 
%% load the position
cfg_pos = [];
cfg_pos.convFact =  [23 33]; % conversion factor from pixels on the camera to cm in the box. I guessed at this one based on Med box dimensions in previous papers. 
pos = LoadPos(cfg_pos); % load the position data

% get the movement speed and periods of immobility
linspeed = getLinSpd([],pos); % linear speed

% smooth the speed a bit
pos_Fs = round(1/mode(diff(linspeed.tvec))); % get the sampling frequency for the speed data. 
linspeed.data = smoothdata(linspeed.data, 'gaussian', pos_Fs * 2 ); % smooth over 2s window with gaussian. 

% Threshold speed. Good for tracking when the are moving but not that
% helpful for LFP analyses since it breaks the data up. 
cfg_l = []; 
cfg_l.method = 'raw'; 
cfg_l.operation = '>'; 
cfg_l.threshold = 5; % speed limit in pixels/sec.  Need to convert to cm/sec to be interpretable. 
iv_fast = TSDtoIV(cfg_l,linspeed); % only keep intervals with speed above thresh

mov_idx = linspeed.data > cfg_l.threshold; % get the data points where the animal is moving. 

% get the nearest pos values for the csc artifacts. 
art_move_idx = unique(nearest_idx3( csc.tvec(art_idx),pos.tvec)); 

%% plot the position

figure(1)
clf % clear the figure

subplot(6,2,[1 3 5])
hold on
plot(pos.data(1,:), pos.data(2,:), '.k')
plot(pos.data(1,art_move_idx), pos.data(2,art_move_idx), '.r')

legend({'position', 'large amp evts'})
title('Tracking')

% or plot it with the time as a colour
subplot(6,2,[2 4 6])
c_ord = winter(length(pos.tvec)); % give it a colour range starting at one color and ending at another. 
scatter(pos.data(1,mov_idx), pos.data(2,mov_idx),3,  c_ord(mov_idx,:)); 
% legend({'position'})
c = colorbar('Location', 'north');
c.Ticks = [0 1];
c.TickLabels = {'start', 'end'};
title('Movement only')

% plot the positions in time when the animal is moving. 
mov_h = subplot(6,2,7:8)
hold on
tvec_zero = pos.tvec - pos.tvec(1);  % makes time realtive to start of recording. 
yyaxis left

plot(tvec_zero, pos.data(1,:), '.k');
plot(tvec_zero, pos.data(2,:), '.k');
% overlay the movement periods in colour
plot(tvec_zero(mov_idx), pos.data(1,mov_idx), '.', 'color', c_ord(1,:));
plot(tvec_zero(mov_idx), pos.data(2,mov_idx), '.', 'color', c_ord(end,:));
xlim([tvec_zero(1), tvec_zero(end)]); % cut the axis limits to the data. Matlab doesn't do this by default anymore
y_lim = ylim; % get the y axis limits

ylim([y_lim(1)-20 max(y_lim(2))])
legend({'x', 'y','x move' 'y move'}, 'Location', 'northoutside', 'Orientation', 'horizontal')
xlabel('time (s)')
ylabel('position (cm)')

yyaxis right
plot(tvec_zero, linspeed.data)
y_lim = ylim; % get the y axis limits

ylim([y_lim(1) max(y_lim(2))+20])
ylabel('speed (cm/s)')

%% Get the theta power in time

% filter the LFP in the theta band
cfg_filt_t = [];
cfg_filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [5 11]; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 3; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool

theta_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using

theta_amp = abs(hilbert(theta_csc.data)); % get the amplitude

theta_phi  = angle(hilbert(theta_csc.data(1,:))); 

theta_csc.data = theta_csc.data(1,:); 
%% Get the gamma power in time

cfg_filt_t = [];
cfg_filt_t.type = 'butter';%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = [30 55]; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 3; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool

gamma_csc = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using

gamma_amp = abs(hilbert(gamma_csc.data)); % get the amplitude

gamma_csc.data = gamma_csc.data(1,:); 
%% Mod idx

mod_th_g = MS_ModIdx_win(theta_csc, gamma_csc, 30*theta_csc.cfg.hdr{1}.SamplingFrequency);

%% plot the theta amp on top of the raw LFP
c_ord = winter(length(csc.label)); 

move_h = subplot(6, 2, 9:12)
hold on
csc_tvec_zero = csc.tvec - csc.tvec(1); % will differ from tvec_zero since this has a higher sampling freq

for ii = length(csc.label):-1:1
plot(csc_tvec_zero, csc.data(ii,:)-(ii/50), 'color', c_ord(ii,:));
plot(csc_tvec_zero(art_idx), csc.data(ii,art_idx)-(ii/50), '.r');

% labels{ii} = csc.label{ii}; 
end

xlim([csc_tvec_zero(1) csc_tvec_zero(end)]); 
xlabel('time (s)')
ylabel('voltage (mV)')
legend(csc.label)
% plot(csc_tvec_zero, theta_amp(1,:), 'b');


% linkaxes(




