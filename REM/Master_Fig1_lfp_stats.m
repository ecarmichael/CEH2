% % Master_Fig1_lfp_stats

lfp_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\LFP';
save_dir =  'C:\Users\ecarm\Desktop\REM_figs';

cd(lfp_dir)
%%  loop over the concatenated sleep LFP data for each session and

l_list = dir('lfp*.h5');

pre_rem_td = []; pre_rem_td_f = [];
pre_sws_td = []; pre_sws_td_f = [];
post_rem_td = []; post_rem_td_f = [];
post_sws_td = []; post_sws_td_f = [];
for ii = length(l_list):-1:1

    if contains(l_list(ii).name, 'pv1254')
        rm_idx(ii) = true;
        continue
    end

    rm_idx(ii) = false;

    lfp = MS_h5_to_stuct(l_list(ii).name);

    % get the theta delta ratio using an overall bandpower across the
    % pre_rem_td(ii) = bandpower(lfp.pre_rem_lfp(2,:), 2000, [6 10])./ bandpower(lfp.pre_rem_lfp(2,:), 2000, [1 4]);
    % pre_sws_td(ii) = bandpower(lfp.pre_sws_lfp(2,:), 2000, [6 10])./ bandpower(lfp.pre_sws_lfp(2,:), 2000, [1 4]);
    % post_rem_td(ii) = bandpower(lfp.post_rem_lfp(2,:), 2000, [6 10])./ bandpower(lfp.post_rem_lfp(2,:), 2000, [1 4]);
    % post_sws_td(ii) = bandpower(lfp.post_sws_lfp(2,:), 2000, [6 10])./ bandpower(lfp.post_sws_lfp(2,:), 2000, [1 4]);

    % continuous ratio version
    % loop for clarity
    for ff = 1:4
            % empty the collectors

        this_td = []; this_td_f = [];

        if ff == 1
            this_data = lfp.pre_rem_lfp(2,:);
            [~, idx] = findpeaks(diff(lfp.pre_rem_lfp(1,:)), 'MinPeakHeight', 20);
        elseif ff == 2
            this_data = lfp.pre_sws_lfp(2,:);
            [~, idx] = findpeaks(diff(lfp.pre_sws_lfp(1,:)), 'MinPeakHeight', 20);
        elseif ff == 3
            this_data = lfp.post_rem_lfp(2,:);
            [~, idx] = findpeaks(diff(lfp.post_rem_lfp(1,:)), 'MinPeakHeight', 20);
        elseif ff == 4
            this_data = lfp.post_sws_lfp(2,:);
            [~, idx] = findpeaks(diff(lfp.post_sws_lfp(1,:)), 'MinPeakHeight', 20);
        end

        % filtered continous ratio
        % theta filter;
        d = fdesign.bandpass('N,F3dB1,F3dB2',4,6,12,2000);
        Hd = design(d,'butter');
        b = Hd.sosMatrix; a = Hd.scaleValues;

        theta = filtfilt(b,a,this_data);
        theta_amp = smoothdata(abs(hilbert(theta)), 'movmean', 2000*10);

        % delta filter
        d = fdesign.bandpass('N,F3dB1,F3dB2',4,1.5,4,2000);
        Hd = design(d,'butter');
        b = Hd.sosMatrix; a = Hd.scaleValues;

        delta = filtfilt(b,a,this_data);
        delta_amp = smoothdata(abs(hilbert(delta)), 'movmean', 2000*10);

        % loop over epcohs to get the bandpower and the median continuous
        % theta delat ratio. 

        if length(idx) == 1
            this_td(end+1) = bandpower(this_data(1,1:idx(1)), 2000, [6 10])./ bandpower(this_data(1,1:idx(1)), 2000, [1 4]);
            this_td_f(end+1) = median(theta_amp(1:idx(1))./delta_amp(1:idx(1)));
        else
            for ib = 1:length(idx)
                if ib == 1
                    this_td(end+1) = bandpower(this_data(1,1:idx(ib)), 2000, [6 10])./ bandpower(this_data(1,1:idx(ib)), 2000, [1 4]);
                    this_td_f(end+1) = median(theta_amp(1:idx(ib))./delta_amp(1:idx(ib))); 

                elseif ib == length(idx)
                    this_td(end+1) = bandpower(this_data(1,idx(ib):end), 2000, [6 10])./ bandpower(this_data(1,idx(ib):end), 2000, [1 4]);
                    this_td_f(end+1) = median(theta_amp(idx(ib):end)./delta_amp(idx(ib):end));
                else
                    this_td(end+1) = bandpower(this_data(1,idx(ib):idx(ib+1)), 2000, [6 10])./ bandpower(this_data(1,idx(ib):idx(ib+1)), 2000, [1 4]);
                    this_td_f(end+1) = median(theta_amp(idx(ib):idx(ib+1))./delta_amp(idx(ib):idx(ib+1)));
                end
                % disp(this_td)
            end
        end

        if ff == 1
            % pre_rem_td_f(ii) = median(theta_amp./delta_amp);
            pre_rem_td_f = [pre_rem_td_f this_td_f]; 
            pre_rem_td = [pre_rem_td this_td];
        elseif ff == 2
            % pre_sws_td_f(ii) = median(theta_amp./delta_amp);
            pre_sws_td_f = [pre_sws_td_f this_td_f];
            pre_sws_td = [pre_sws_td this_td];
        elseif ff == 3
            % post_rem_td_f(ii) = median(theta_amp./delta_amp);
            post_rem_td_f = [post_rem_td_f this_td_f];
            post_rem_td = [post_rem_td this_td];
        elseif ff == 4
            % post_sws_td_f(ii) = median(theta_amp./delta_amp);
            post_sws_td_f = [post_sws_td_f this_td_f];
            post_sws_td = [post_sws_td this_td];
        end
    end

    clearvars lfp
end

%%
sws_c = [132 147 35]/255;
rem_c = [67 127 151]/255;

figure(1001)
clf

% using the continuous filtered signals.
subplot(2,2,1)
[~, ~, ~, p, stats] = MS_bar_w_err([pre_rem_td_f post_rem_td_f], [pre_sws_td_f, post_sws_td_f], [rem_c; sws_c], 1, 'ttest2', 1:2);
ylabel({'theta/delta'; 'continuous'})
set(gca, 'XTickLabel', {'REM' 'NREM'})
xlim([0 3])
axis square
fprintf('REM (%0.2f +/- %0.2f ) - NREM (%0.2f +/- %0.2f ) theta delta ration: t(%d) = %0.2f; p %0.3f\n',...
    mean([pre_rem_td_f post_rem_td_f]), MS_SEM([pre_rem_td_f post_rem_td_f]),...
    mean([pre_sws_td_f, post_sws_td_f]), MS_SEM([pre_sws_td_f, post_sws_td_f]),...
    stats.df, stats.tstat, p )


% using bandpower per epoch 
subplot(2,2,3)
[~, ~, ~, p, stats] = MS_bar_w_err([pre_rem_td post_rem_td], [pre_sws_td, post_sws_td], [rem_c; sws_c], 1, 'ttest2', 1:2);
ylabel({'theta/delta'; 'bandpower'})
set(gca, 'XTickLabel', {'REM' 'NREM'})
axis square
xlim([0 3])
fprintf('REM (%0.2f +/- %0.2f ) - NREM (%0.2f +/- %0.2f ) theta delta ration: t(%d) = %0.2f; p %0.3f\n',...
    mean([pre_rem_td post_rem_td]), MS_SEM([pre_rem_td post_rem_td]),...
    mean([pre_sws_td, post_sws_td]), MS_SEM([pre_sws_td, post_sws_td]),...
    stats.df, stats.tstat, p )

% set figure properties
cfg_fig.ft_size = 7; 
set(gcf, 'Position', [.3 .1 .3 .4])
SetFigure(cfg_fig, gcf); 

%% save
print("-bestfit",[save_dir filesep 'fig1_theta_delta'], '-dpdf', "-vector")
