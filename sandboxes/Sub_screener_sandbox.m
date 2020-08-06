
%% Sub screening sandbox

close all
restoredefaultpath
global PARAMS  % these are global parameters that can be called into any function.  I limit these to directories for storing, loading, and saving files and codebases.
os = computer;

if ismac
    %     PARAMS.data_dir = '/Users/jericcarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    %     PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    %     PARAMS.stats_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    %     PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    %     PARAMS.code_CEH2_dir = '/Users/jericcarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    %
elseif strcmp(os, 'GLNXA64')
    
    %     %     PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    %     PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/7_12_2019_PV1069_LTD5'; % where to find the raw data
    %     %     PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/EV/';
    %     PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/'; % raw data location.
    %     PARAMS.inter_dir = '/home/ecarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    %     PARAMS.stats_dir = '/home/ecarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    %     PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    %     PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    %
else
    PARAMS.data_dir = 'J:\Williams_Lab\II_classification'; % where to find the raw data
    PARAMS.raw_data_dir = 'J:\Williams_Lab\II_classification'; % raw data location.
    PARAMS.inter_dir = 'J:\Williams_Lab\II_classification\Inter'; % where to put intermediate files
    PARAMS.stats_dir = 'J:\Williams_Lab\II_classification\Inter\Stats'; % where to put the statistical output .txt
    PARAMS.code_base_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2'; % where the multisite repo can be found
    PARAMS.code_seqnmf_dir = 'C:\Users\ecarm\Documents\GitHub\seqNMF'; % where the multisite repo can be found
    
end

% colours
PARAMS.L_grey = [0.8 0.8 0.8];
PARAMS.D_grey = [0.2 0.2 0.2];
PARAMS.blue = [0.3639    0.5755    0.7484];
PARAMS.red = [0.9153    0.2816    0.2878];
PARAMS.green= [0.4416    0.7490    0.4322];
PARAMS.gold = [1.0000    0.5984    0.2000];

rng(11,'twister') % for reproducibility


% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_CEH2_dir));
cd(PARAMS.raw_data_dir) % move to the data folder

% make the Inter and Stats dir if they don't exist.
if ~exist(PARAMS.stats_dir,'dir')
    mkdir(PARAMS.stats_dir)
end

% % try the newer NLX loaders for UNIX
% [~, d] = version;
% if str2double(d(end-3:end)) >2014 && strcmp(os, 'GLNXA64')
%     rmpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared/io/neuralynx')
%     addpath(genpath('/Users/jericcarmichael/Documents/NLX_loaders_UNIX_2015'))
%     disp('Version is greater than 2014b on UNIX so use updated loaders found here:')
%     which Nlx2MatCSC
% end

clear d os

%%  load data

load('567_ms.mat')
load('567_behav.mat')
ms = MS_msExtractBinary_detrendTraces(ms, 2);

%% try to classify the cells as long or short tailed
peak_window = 33; % how far from the onset of the transient do you look for the peak.
tail_thresh = 3.5; % try to match Roy et al using 3.5s.  Need to convert to samples using dT.
tail_thresh_idx = ceil(tail_thresh/mode(diff(ms.time))); % Need to convert to samples using dT.
% sample_window = tail_thresh*4; % number of samples from the transient onset to look for the high_thresh. Used 4*tail_thresh since anything longer than tail_thresh will be classified as long any ways.
high_thresh = 0.33; % percent of local maxima as cut off for that transient. Again from Roy 2017.
all_trans_dur = []; % array for all the durations across cells/transients. used for distribution plot a la Roy 2017 Fig S6c
min_separation = 0.250; % minimum distance from a previous/subsequent event.  From Roy 2017;


frame_time = 0.001*(mode(diff((ms.time))));  % convert time frame to seconds.
cfg.plot = 0;
cfg.overlap = 0;

for iC = size(ms.RawTraces,2):-1:1  % loop over cells in the ms struct
    fprintf('Processing cell #%d...', iC);
    these_dur = [];
    if sum(ms.Binary(:,iC)) ==0  % check for any transients passing threshold.  If none, continue.
        all_dur(iC) = NaN;
        continue
    else
        % Get the transient onset idx
        [pks, loc] = findpeaks(diff(ms.Binary(:,iC)), 'MinPeakDistance',33);
        
        overlap = zeros(1,length(loc)); % hold a record of transients that overlapped with prior transients.
        
        diff_loc_sec = diff(loc)*frame_time;
        loc =loc+1; %offset diff
        if cfg.plot
            
            figure(11) % plot for debugging
            hold on
            h(1) =  plot(ms.time, ms.detrendRaw(:,iC));
            h(2) =  plot(ms.time, ms.Binary(:,iC));
            h(3) =  plot(ms.time(loc), pks, '*');
            h(4) =  plot(ms.time(2:end), diff(ms.Binary(:,iC)));
        end
        %
        %
        for iTran = 1:length(pks) % loop transients in this Ca trace
            
            
            % find the max within a specified window
            if loc(iTran)+peak_window <= length(ms.detrendRaw(:,iC));
                [val, idx] = max(ms.detrendRaw(loc(iTran):loc(iTran)+peak_window,iC)); % get the peak val of the maxima.  Use the idx to start search for
            else
                [val, idx] = max(ms.detrendRaw(loc(iTran):end,iC)); % if beyond length of recording, then just use the end. this will be flagged as an overlap.
            end
            
            idx = idx+loc(iTran);
            
            %             off_val = val-ms.RawTraces(loc(iTran)); % correct for starting value of the transient
            %             val-(off_val*high_thresh);
            if iTran == length(pks) % check if this is the last one in the set. If so use last sample.
                
                end_idx = find(ms.detrendRaw(idx:end,iC) <= val*high_thresh);
                
                if cfg.plot
                    hold on
                    h(5) =plot(ms.time(loc(iTran):loc(iTran)+peak_window),ms.detrendRaw(loc(iTran):loc(iTran)+peak_window,iC), 'm');
                    if ~isempty(end_idx)
                        h(6) =plot(ms.time(idx:idx+end_idx(1)),ms.detrendRaw(idx:idx+end_idx(1),iC), '-c');
                    else
                        h(6) =  plot(ms.time(idx:end),ms.detrendRaw(idx:end,iC), '-c');
                    end
                    %                 xlim([ms.time(loc(iTran)), ms.time(loc(iTran+1))])
                    h(7) = hline(val*high_thresh);
                end
                if ~cfg.overlap && isempty(end_idx) % toggle if the transients are allowed to overlap with the next or the end.
                    end_idx = NaN;
                    overlap(iTran) = 1;
                elseif cfg.overlap && isempty(end_idx)
                    end_idx = length(ms.detrendRaw) -loc(iTran); % unlikely event that it does not end before the end.
                    overlap(iTran) = 1;
                end
            else
                end_idx = find(ms.detrendRaw(idx:loc(iTran+1),iC) <= val*high_thresh);
                if cfg.plot
                    hold on
                    h(5) = plot(ms.time(loc(iTran):loc(iTran)+peak_window),ms.detrendRaw(loc(iTran):loc(iTran)+peak_window,iC), 'm');
                    if ~isempty(end_idx)
                        h(6) = plot(ms.time(idx:idx+end_idx),ms.detrendRaw(idx:idx+end_idx,iC), '-c');
                    else
                        h(6) = plot(ms.time(idx:loc(iTran+1)),ms.detrendRaw(idx:loc(iTran+1),iC), '-c');
                    end
                    %                 xlim([ms.time(loc(iTran)), ms.time(loc(iTran+1))])
                    h(7) = hline(val*high_thresh);
                end
                if ~cfg.overlap && isempty(end_idx)
                    end_idx = NaN;
                    overlap(iTran) = 1;
                elseif cfg.overlap && isempty(end_idx)
                    end_idx = loc(iTran+1) -loc(iTran);
                    overlap(iTran) = 1;
                end
            end
            
            if iTran > 1 && overlap(iTran-1) % check if the previous transient didn't reach the threshold before this transient began.
                end_idx = NaN;
            end
            %             end_idx = end_idx(1); % get the first index where it drops below the threshold.
            
            
            
            these_dur(iTran) = end_idx(1)*frame_time; % get the number of samples before it crosses the threshold.  Convert to time w/ dT. Gives time in Seconds
            
        end % end transient loop.
        
        
        all_dur(iC) = nanmean(these_dur); % get the average duration for across transients for this cell.
        
        all_trans_dur = [all_trans_dur, these_dur]; % collect all the durations for plotting a distribution a la
    end % above binary threshold check.
    fprintf('done.\n')
    if cfg.plot;  cla(h); end % clear the plot axes.
end


%% make a plot of the distributions
n_bins = 50;

figure(101)
log_dur = log10(all_dur); % convert to log10
log_thres = log10(tail_thresh);  % match Roy 2017

[N, X] = hist(log_dur, n_bins);
bar(X, N/length(all_dur));
vline(log_thres)
xlabel('log Ca event duration (s)')
ylabel('frequency')
title('Mean Ca^2^+ transient duration across all cells')
text(log_thres,max(N/length(all_dur))*0.9,[num2str(tail_thresh) 's cutoff'],'HorizontalAlignment', 'left')
xlim([-1 4]);
[mval, midx] = max(N/length(all_dur)); % find the peak idx for text.
text(X(midx),mval,[num2str(10^X(midx), '%2.2f') 's mean'],'HorizontalAlignment', 'left')


axes('Position',[.7 .7 .2 .2])
box on
[N, X] = hist(all_dur, n_bins);
bar(X, N/length(all_dur));
vline(tail_thresh);
xlabel('Raw Ca^2^+ duration (s)')

xlim([0 20]);

%% align behaviour and Ca

if behav.time(end) ~= ms.time(end) || length(behav.time) ~= length(ms.time)
    fprintf('<strong> %s </strong>: behaviour and Ca are not the same length ro end time.  attempting alignment \n', mfilename);
    
    
    behav_aligned = MS_align_data(behav, ms);
%     behav_aligned = MS_align_data(behav, ms);

end

left_idx = MS_get_direction(behav_aligned.position(:,1), -.1); % use -threshold for leftbound and + for right.
right_idx = MS_get_direction(behav_aligned.position(:,1), .1); % use -threshold for leftbound and + for right.


movement_thresh = 2.5; % in cm/s
movement_idx = behav_aligned.speed >movement_thresh; % get times when the animal was moving.
left_idx = left_idx & movement_idx; % only keep the indices when they are moving to the left.
right_idx = right_idx & movement_idx;

[L_laps, L_lap_start_idx, L_lap_end_idx] = MS_get_laps(left_idx, floor(1.5*(1/frame_time)),floor(10*(1/frame_time)));
[R_laps, R_lap_start_idx, R_lap_end_idx] = MS_get_laps(right_idx, floor(1.5*(1/frame_time)),floor(10*(1/frame_time)));

%% plot basics
for iC = 1; % loop through cells
    
    M = 7; % rows
    N = 6; % columns
    figure(100) % main figure
    
    %%% title information
    subplot(M, N, [N-2 (N*2)-2]) % title information. Top right corner,
    text(0, 10, ['Cell id: ' num2str(iC)])
    text(0, 8, ['Binary thresh: ' num2str(ms.Binary_threshold)])
    
    if all_dur(iC) >= tail_thresh
        text(0, 6, ['Mean tail length: ' num2str(all_dur(iC),'%2.2f') 's "Long"'])
    else
        text(0, 6, ['Mean tail length: ' num2str(all_dur(iC),'%2.2f') 's "Short"'])
    end
    ylim([0 10])
    axis off
    
    
    %%% raw trace
    subplot(M, N, 1:3)
    plot(ms.time/1000, ms.RawTraces(:,iC), 'color', PARAMS.blue)
    xlim([ms.time(1)/1000 ms.time(end)/1000]);
    xlabel('time(s)');
    
    
    
    %%% X Y position
    subplot(M, N, N+1:N+3)
    hold on
    plot(behav_aligned.time/1000, behav_aligned.position(:,1), 'color', PARAMS.L_grey)
    plot(behav_aligned.time(left_idx)/1000, behav_aligned.position(left_idx,1),'.', 'color', PARAMS.green)
    plot(behav_aligned.time(right_idx)/1000, behav_aligned.position(right_idx,1),'.', 'color', PARAMS.blue)
%     plot(behav_aligned.time/1000, L_laps*10,'.-', 'color', PARAMS.gold, 'linewidth', 0.25)
%     plot(behav_aligned.time/1000, R_laps*10,'.-', 'color', PARAMS.gold,'linewidth', 0.25)

    %
    % plot(behav_aligned.time(lap_start_idx)/1000, behav_aligned.position(lap_start_idx,1),'d', 'color', PARAMS.green)
    % plot(behav_aligned.time(lap_end_idx)/1000, behav_aligned.position(lap_end_idx,1),'d', 'color', PARAMS.green)
    
    
    % plot(behav_aligned.time/1000, behav_aligned.position(:,2),'color', PARAMS.blue)
    ylabel('delta position')
    xlim([behav_aligned.time(1)/1000 max(behav_aligned.time)/1000]);
    % legend({'x', 'y'})
    
    
    %%% speed info
    subplot(M, N, N*2+1:N*2+3)
    hold on
    plot(behav_aligned.time/1000, behav_aligned.speed, 'color', PARAMS.L_grey)
    plot(behav_aligned.time(movement_idx)/1000, behav_aligned.speed(movement_idx),'.', 'color', PARAMS.gold, 'markersize', 1)
    % legend('Speed', 'box', 'off')
    
    xlim([behav_aligned.time(1)/1000 max(behav_aligned.time)/1000]);
    ylabel('Speed cm/s')
    xlabel('time (s)')
    
    
    % %%% orientation info
    % subplot(M, N, N*3+1:N*3+3)
    % plot(behav_aligned.time/1000,ones(size(behav_aligned.time)), 'color', 'w')
    % hold on
    % text(behav_aligned.time(floor(length(behav_aligned.time)/3))/1000, pi, 'HD placeholder')
    % ylabel('HD')
    % ylim([-pi pi])
    % set(gca, 'ytick', [-pi pi], 'yticklabel', {'-pi' 'pi'})
    
    %%% plot the binary times on the position
    subplot(M, N, [ N*3+4:N*3+6]) % N*4+4:N*4+6
    hold on
    plot(behav_aligned.position(:,1), behav_aligned.position(:,2), 'color', PARAMS.L_grey)
    xlim([min(behav_aligned.position(:,1)) max(behav_aligned.position(:,1))])
    ylim(round([min(behav_aligned.position(:,2)) max(behav_aligned.position(:,2))]));
    %get binary 'event times'
    t_binary = find(ms.Binary(:,iC) ==1);
    % put dots on positions when the cell was active.
    
    plot(behav_aligned.position(t_binary,1), behav_aligned.position(t_binary,2),'.', 'color', PARAMS.red)
    xlabel('position (cm)');
    ylabel('position (cm)');
    set(gca, 'ytick', round([min(behav_aligned.position(:,2)) max(behav_aligned.position(:,2))]));
    
    
    % get the transient/position values
    % tran_x = interp1(behav_aligned.time(1:end-1),behav_aligned.position(1:end-1,1),ms.time(t_binary),'linear');
    % tran_y = interp1(behav_aligned.time(1:end-1),behav_aligned.position(1:end-1,2),ms.time(t_binary),'linear');
    %
    % plot(tran_x,tran_y,'.', 'color', PARAMS.red);
    
    
    %%% update position in time with binary 'spikes'
    subplot(M, N, N+1:N+3)
    plot(behav_aligned.time(t_binary)/1000,behav_aligned.position(t_binary,1),'.', 'color', PARAMS.red);
    plot(behav_aligned.time(t_binary)/1000,behav_aligned.position(t_binary,2),'.', 'color', PARAMS.red);
    
    %%% update speed in time with binary 'spikes'
    subplot(M, N, N*2+1:N*2+3)
    plot(behav_aligned.time(t_binary)/1000,behav_aligned.speed(t_binary,1),'.', 'color', PARAMS.red);
    
    
    %%% add the SPF for this cell.
    subplot(M, N, [N (N*2)]) % spf with centroid.
    imagesc(ms.SFPs(:,:,iC))
    hold on
    
    
    %%% Plot the laps
    subplot(M, N, N*5+4:N*5+6)
    MS_plot_laps(behav_aligned, R_laps, t_binary)
    ylabel('R laps');
    
    subplot(M, N, N*6+4:N*6+6)
    MS_plot_laps(behav_aligned, L_laps, t_binary)
    ylabel('L laps');
    
    %%% add in the place/spatial information?
    subplot(M, N, N*4+4:N*4+6)
    bins = min(behav_aligned.position(:,1)):2.5:max(behav_aligned.position(:,1));
    bin_centers = bins + 2.5/2;
    bin_centers = bin_centers(1:end-1);
    [MI, posterior, occupancy, p_active, likelihood] = MS_get_spatial_information(ms.Binary(:,iC), ms.time, behav_aligned.position(:,1), bins);
    plot(bin_centers, likelihood)
    ylabel('p active')
    
    %%% add speed mod score
    
    
    %%% customize figure stuff
    
    pos = get(gcf, 'position');
    set(gcf, 'position', [pos(1)-pos(1)*.6 pos(2)-pos(2)*.6 pos(3)*1.8 pos(4) *1.6])
    
    
    pause(1)
    close(100)
end







