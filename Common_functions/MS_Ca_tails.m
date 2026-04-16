function ms = MS_Ca_tails(ms, t_samples)




%% initialize

cfg.plot =0;
cfg.overlap = 0; 

peak_window = 33; % how far from the onset of the transient do you look for the peak.
tail_thresh = 3.5; % try to match Roy et al using 3.5s.  Need to convert to samples using dT.
tail_thresh_idx = ceil(tail_thresh/t_samples); %ceil(tail_thresh/mode(diff(ms.time))); % Need to convert to samples using dT.
% sample_window = tail_thresh*4; % number of samples from the transient onset to look for the high_thresh. Used 4*tail_thresh since anything longer than tail_thresh will be classified as long any ways.
high_thresh = 0.33; % percent of local maxima as cut off for that transient. Again from Roy 2017.
all_trans_dur = []; % array for all the durations across cells/transients. used for distribution plot a la Roy 2017 Fig S6c
min_separation = 0.250; % minimum distance from a previous/subsequent event.  From Roy 2017;

tvec = 0:length(ms.RawTraces)-1;

%% loop over cells

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

        % diff_loc_sec = diff(loc)*frame_time;
        loc =loc+1; %offset diff
        if cfg.plot

            figure(11) % plot for debugging
            clf
            hold on
            h(1) =  plot(tvec, ms.FiltTraces(:,iC));
            h(2) =  plot(tvec, ms.Binary(:,iC));
            h(3) =  plot(tvec(loc), pks, '*');
            h(4) =  plot(tvec(2:end), diff(ms.Binary(:,iC)));
        end

        for iTran = 1:length(pks) % loop transients in this Ca trace

            % find the max within a specified window
            if loc(iTran)+peak_window <= length(ms.detrendRaw(:,iC))
                [val, idx] = max(ms.detrendRaw(loc(iTran):loc(iTran)+peak_window,iC)); % get the peak val of the maxima.  Use the idx to start search for
            else
                [val, idx] = max(ms.detrendRaw(loc(iTran):end,iC)); % if beyond length of recording, then just use the end. this will be flagged as an overlap.
            end
            
            idx = idx+loc(iTran);

            if idx >= length(ms.detrendRaw)
                continue
            end

            %             off_val = val-ms.RawTraces(loc(iTran)); % correct for starting value of the transient
            %             val-(off_val*high_thresh);
            if iTran == length(pks) % check if this is the last one in the set. If so use last sample.

                end_idx = find(ms.detrendRaw(idx:end,iC) <= val*high_thresh);

                if cfg.plot
                    hold on
                    h(5) =plot(tvec(loc(iTran):end),ms.detrendRaw(loc(iTran):end,iC), 'm');

                        h(6) =  plot(tvec(idx:end),ms.detrendRaw(idx:end,iC), '-c');
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
                    h(5) = plot(tvec(loc(iTran):loc(iTran)+peak_window),ms.detrendRaw(loc(iTran):loc(iTran)+peak_window,iC), 'm');
                    if ~isempty(end_idx)
                        h(6) = plot(tvec(idx:idx+end_idx(1)),ms.detrendRaw(idx:idx+end_idx(1),iC), '-c');
                    else
                        h(6) = plot(tvec(idx:loc(iTran+1)),ms.detrendRaw(idx:loc(iTran+1),iC), '-c');
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



            these_dur(iTran) = end_idx(1)*t_samples; % get the number of samples before it crosses the threshold.  Convert to time w/ dT. Gives time in Seconds

        end % end transient loop.


        all_dur(iC) = median(these_dur, 'omitmissing'); % get the average duration for across transients for this cell.

        all_trans_dur = [all_trans_dur, these_dur]; % collect all the durations for plotting a distribution a la
    end % above binary threshold check.
    fprintf('done.\n')
    if cfg.plot;  cla(h); end % clear the plot axes.
end

%% make a plot of the distributions
n_bins = 50;
c_ord = MS_linspecer(3); 
% 
figure(102)
% subplot(3,1,1:2)
% cfg_plot.view =[0 75];
% cfg_plot.plot_type = '2d';
% MS_plot_ca(cfg_plot, ms)



% subplot(3,1,3)
log_dur = log10(all_dur); % convert to log10
log_thres = log10(tail_thresh);  % match Roy 2017

[N, X] = hist(log_dur, -1:.025:log10(5));
bar(X, N/length(all_dur), 'FaceColor', c_ord(1,:), 'EdgeColor','none');
xline(log_thres)
% xline(log10(.25), 'k', '0.25s', 'FontSize',14)
% xline(log10(.5),'k', '0.5s', 'FontSize',14)
% xline(log10(1), 'k', '1s', 'FontSize',14)
xlabel('Ca event duration (s)')
ylabel('frequency')
title('Mean Ca^{2+} transient duration across all cells','Interpreter', 'Tex')
text(log_thres,max(N/length(all_dur))*1,[num2str(tail_thresh) 's cutoff'],'HorizontalAlignment', 'left', 'FontSize',14)
xlim([-1 log10(5)]);
[mval, midx] = max(N/length(all_dur)); % find the peak idx for text.
text(X(midx),mval*1.1,[num2str(10^X(midx), '%2.2f') 's median'],'HorizontalAlignment', 'left', 'FontSize',14)
ylim([0 .15])
set(gca, "Xtick", log10([.1 .25 .5 .75 1 2.5 5 10]), 'XTickLabel',[.1 .25 .5 .75 1 2.5 5 10]);
SetFigure([], gcf, 1)
