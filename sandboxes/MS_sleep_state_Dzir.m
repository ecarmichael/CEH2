function MS_sleep_state(cfg_in, csc)
%% MS_sleep_state: WIP sleep state detector.  Currently based on Dzirasa et al. 2006
%
%
%
%    Inputs:
%     - csc [struct] in the MS_loadCSC format.
%
%
%
%    Outputs:
%     - none ATM
%
%
%
%
% EC 2020-02-17   initial version
%
%
%% Initialize
cfg_def = [];
cfg_def.f_bands = [2,2,2,2; 4.5, 9 ,20,55]'; % which frequencies to look for
cfg_def.lfp_chan = 2; % which channel in the csc to use for the lfp
cfg_def.emg_chan = 1; % which channel in the csc to use for the emg.
cfg_def.emg_f_range = 100:2:200; % used to compute spectrogram.
cfg_def.lfp_f_range = 1:.5:64; % similar to Dzirasa et al. 2006
cfg_def.ratios = [1,3,4;2,4,5]'; % frequency ratios to use.  default is 2-4.5 ./ 2-9, 2-20 ./2-55, 2-55 ./ emg
cfg_def.gau_win = 30;% how many samples for smoothing.  since the spectrogram computes over 1s windows this values is in seconds.
cfg_def.cut = []; % use if you want to cut out data.
cfg_def.norm =0; % normalize the power to the max value for the LFP and the EMG. default is 0 (no)

cfg = ProcessConfig(cfg_def, cfg_in);

%% get some sleep ratios and clustering? based on Dzirasa et al. 2006 (mice HC sleep states)
% get the ratios and then bin the data to get discrete points. maybe 1s or
% 5sec?  Then run kmeans on those clusters to pull out W, SW, and REM.
% Maybe add in a statement that says if in SW and there is a short W
% followed by long SW call this a micro Arrosal .


%% use spectrogram approach. (not great...)
figure(101)
plot(csc.tvec, csc.data(cfg.lfp_chan,:))
title(['LFP Channel: ' csc.label(cfg.lfp_chan) ' ' num2str(cfg.lfp_chan)])
if cfg.cut ==1
    fprintf('Select data to remove.  only one block ATM...')
    [cuts, ~] = ginput(2);
    rec_1 = restrict(csc, csc.tvec(1), cuts(1));
    rec_2 = restrict(csc, cuts(2), csc.tvec(end));
    close(101)
    temp_data = [rec_1.data, rec_2.data];
else
    temp_data = csc.data;
end

% temp_tvec = 0:1/csc.cfg.hdr{1}.SamplingFrequency:length(temp_data)/csc.cfg.hdr{1}.SamplingFrequency;
% temp_tvec = temp_tvec(1:end-1);
figure(1111)
% get the LFP spectrogram
[~,F,T,P] = spectrogram(temp_data(cfg.lfp_chan,:), hanning(csc.cfg.hdr{2}.SamplingFrequency*2), csc.cfg.hdr{2}.SamplingFrequency, cfg.lfp_f_range, csc.cfg.hdr{2}.SamplingFrequency);
subplot(2,1,1)
imagesc(T,F,10*log10(P));
axis xy; ylabel('freq'); xlabel('time(s)');
title(csc.label{cfg.lfp_chan})

% get the emg spectrogram.
[~, F_emg, T_emg, P_emg] = spectrogram(temp_data(cfg.emg_chan,:), hanning(csc.cfg.hdr{2}.SamplingFrequency*2), csc.cfg.hdr{2}.SamplingFrequency, cfg.emg_f_range, csc.cfg.hdr{1}.SamplingFrequency);
subplot(2,1,2)
imagesc(T_emg,F_emg,10*log10(P_emg));
axis xy; ylabel('freq'); xlabel('time(s)');
title(csc.label{cfg.emg_chan})


% f_label = {'low1', 'low2', 'wide_low', 'wide'};

for iF = 1:length(cfg.f_bands)
    f_label{iF} = ['f_' num2str(cfg.f_bands(iF,1)) '_' strrep(num2str(cfg.f_bands(iF,2)), '.', 'p')];
    
    f_idx = find(cfg.f_bands(iF,1) <= F & F <= cfg.f_bands(iF,2));
    if cfg.norm ==1
        pow.(f_label{iF}) = (mean(10*log10(P(f_idx,:))))./ max(mean(10*log10(P(f_idx,:))));
    else
        pow.(f_label{iF}) = mean(10*log10(P(f_idx,:)));
    end
    
end
if cfg.norm ==1
    pow.emg = (mean(10*log10(P_emg)))./max(mean(10*log10(P_emg)));
else
    pow.emg = mean(10*log10(P_emg));
end

f_label{end+1} = 'emg';

% made it to here.  Don't know if this works.
ratios_cat = [];
for iR = 1:length(cfg.ratios)
    ratios{iR} = pow.(f_label{cfg.ratios(iR,1)}) ./ pow.(f_label{cfg.ratios(iR,2)});
    ratio_label{iR} = [f_label{cfg.ratios(iR,1)} '/' f_label{cfg.ratios(iR,2)}];
    
    ratios_con{iR} = smoothdata(ratios{iR},'gaussian',cfg.gau_win);
    ratios_cat = [ratios_cat, ratios_con{iR}'];
end

% ratio_2 = pow.(f_label{1}) ./ pow.(f_label{2}) ;
% ratio_1 = pow.wide_low ./ pow.wide;
% ratio_t_emg = pow.wide ./ emg;

% ratio_1_con = conv2(ratio_1, gausswin(csc.cfg.hdr{1}.SamplingFrequency*20),'same');
% ratio_2_con = conv2(ratio_2, gausswin(csc.cfg.hdr{1}.SamplingFrequency*20),'same');
% ratio_t_emg_con = conv2(ratio_t_emg, hanning(csc.cfg.hdr{1}.SamplingFrequency*20),'same');
%
% ratio_1_con =  smoothdata(ratio_1,'gaussian',30);
% ratio_2_con =  smoothdata(ratio_2,'gaussian',30);
% ratio_t_emg_con =  smoothdata(ratio_t_emg,'gaussian',30);

% combine for kmeans
% ratios_12 = [ratio_1_con', ratio_2_con', ratio_t_emg_con'];

[idx,C] = kmeans(ratios_cat,3);

figure(1112)
% plot the power
ax_s(1) = subplot(2, 4, 1:2);
hold on
for iF = 1:length(f_label)
    plot(T_emg, pow.(f_label{iF}));
end
legend(strrep(f_label, '_', '-'))
if cfg.norm == 1
    ylabel('norm power')
else
   ylabel('raw power') 
end
% plot the ratios
ax_s(2) = subplot(2, 4, 5:6);
hold on
for iR = 1:length(ratio_label)
    plot(T_emg, ratios_con{iR});
end
ylabel('Ratios')
legend(strrep(ratio_label, '_', '-'))

linkaxes(ax_s, 'x')
subplot(2, 4, [3,4,7,8])
MS_gscatter3d(ratios_cat,idx,linspecer(3));
hold on
plot3(C(:,1),C(:,2),C(:,3),'kx', 'markersize', 12, 'linewidth', 5)
legend('Cluster 1','Cluster 2','Cluster 3','Cluster Centroid');
xlabel(strrep(ratio_label{1}, '_', '-'));
ylabel(strrep(ratio_label{2}, '_', '-'));
zlabel(strrep(ratio_label{3}, '_', '-'));
view(3)

%% filtered LFP method: not working well.
% % delta filter.
% cfg_filt_d = [];
% cfg_filt_d.type = 'fdesign'; %the type of filter I want to use via filterlfp
% cfg_filt_d.f  = [2 4.5];
% cfg_filt_d.order = 8; %type filter order
% cfg_filt_d.display_filter = 1; % use this to see the fvtool and wait for a keypress before continuing.
% low1_csc = FilterLFP(cfg_filt_d,csc);
% close all
%
%
% % filter into the theta band
% cfg_filt_t = [];
% cfg_filt_t.type = 'cheby1';%'fdesign'; %the type of filter I want to use via filterlfp
% cfg_filt_t.f  = [2 9];
% cfg_filt_t.order = 3; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% % reason.  .
% low2_csc = FilterLFP(cfg_filt_t, csc);
%
% % 'wide' from Watson et al. 2016
% % filter into the theta band
% cfg_filt_t = [];
% cfg_filt_t.type = 'fdesign'; %the type of filter I want to use via filterlfp
% cfg_filt_t.f  = [2 20];
% cfg_filt_t.order = 12; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% % reason.  .
% low_wide_csc = FilterLFP(cfg_filt_t, csc);
%
%
% % 'wide' from Watson et al. 2016
% % filter into the theta band
% cfg_filt_t = [];
% cfg_filt_t.type = 'fdesign'; %the type of filter I want to use via filterlfp
% cfg_filt_t.f  = [2 55];
% cfg_filt_t.order = 8; %type filter order
% cfg_filt_t.display_filter = 1; % use this to see the fvtool (but very slow with ord = 3 for some
% % reason.  .
% wide_csc = FilterLFP(cfg_filt_t, csc);
%
% ratio_2 = abs(hilbert(low_wide_csc.data(2,:))) ./ abs(hilbert(wide_csc.data(2,:)));
% ratio_1 = abs(hilbert(low1_csc.data(2,:))) ./ abs(hilbert(low2_csc.data(2,:)));
%
% % ratio_con_2 = conv
%
% % ratio_1_2CSC = ratio

