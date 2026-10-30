function S = MS_get_wave_properties_phy2(S, waves,win_tvec, plot_flag)
%% MS_get_wave_properties_phy2:
%
%
%
%    Inputs:
%    - S [struct]   Spike 'S' structure from LoadSpikes
%
%    - wave  [5 x nSamples]  waveform array.
%
%    - plot_flag: [logical]  0 = no plots, 1 = plots.
%
%    Outputs:
%    - S [struct]  TS structure with the features added to the 'usr' field.
%
%
%
%
%
% EC 2021-07-22   initial version
% EC 2026-10-29:  PHy2 version based on https://www.cell.com/neuron/fulltext/S0896-6273(25)00665-8#sec-8
%
%
%% initialize
if nargin < 3
    plot_flag = 0;
end

%% loop over cells


for ii = length(S.t):-1:1

    log_isi = log10(mean(diff(S.t{ii})));
    ISI = diff(S.t{ii});

    % %% get the autocorrelation (from vandermeerlab.org wiki: https://rcweb.dartmouth.edu/~mvdm/wiki/doku.php?id=analysis:nsb2016:week9)
    % xbin_centers = -.025-0.001:0.001:0.025+0.001; % first and last bins are to be deleted later
    % ac = zeros(size(xbin_centers));
    %
    % for iSpk = 1:length(S.t{1})
    %
    %     relative_spk_t = S.t{1} - S.t{1}(iSpk);
    %
    %     ac = ac + hist(relative_spk_t,xbin_centers); % note that hist() puts all spikes outside the bin centers in the first and last bins! delete later.
    %
    % end
    %
    % xbin = xbin_centers(2:end-1); % remove unwanted bins
    % zero_idx = find(xbin == 0);
    % ac = ac(2:end-1);
    % ac(zero_idx) = 0;
    %
    %
    %
    % % ac = ac./max(ac); % normalize
    % if plot_flag
    %     figure
    %     clf;
    %     subplot(2,3,3)
    %     bar(xbin*1000, ac)
    %     xlabel('Lag (ms)');
    %
    %     subplot(2,3,6)
    %     histogram(log10(ISI),50,'Normalization', 'probability')
    %     % set(gca, 'XScale', 'log')
    %     xlabel('ISI (sec)');
    %     xline(prctile(log10(ISI), 10), '--', '10th prctl')
    %     xlim(log10([0.0001 1]))
    %     set(gca, 'XTickLabel', 10.^get(gca, 'xtick'))
    %
    % end
    %
    % bin_idx_2_6 = zero_idx +2: zero_idx+6;
    % bin_idx_2_20 = zero_idx +2: zero_idx+20;



    % bursting propensity
    burst_idx = sum(ISI < .010) / sum(ISI > .010);


    %% wave features
    wave = waves(ii,:);

    %key time points
    z_idx = nearest_idx(0, win_tvec);
    one_idx =  nearest_idx(1, win_tvec);
    neg_idx = nearest_idx(-.25, win_tvec); 

    % trough
    [t, t_idx] = min(wave);

    % (post-trough) peak
    [p, p_idx] = max(wave(z_idx:one_idx));
    p_idx = p_idx+z_idx;

    % pre peak
    [pre, pre_idx] = max(wave(neg_idx:z_idx));
    pre_idx = pre_idx+(neg_idx-1); 

    % duration
    dur = win_tvec(p_idx) - win_tvec(t_idx);

    % amplitude
    amp = abs(p) - abs(t); 

    % pre- peak symmetry ('ab_ratio' in cell explorer)
     pre_peak_r =  (p-pre)/(pre+p); 

    % trough - peak ratio
    troughPeakRatio = abs(p) / abs(t);

    % repolarization slope
    repolSlope = (p - t) / (win_tvec(p_idx) - win_tvec(t_idx));

    % pre-peak to trough ratio
    prePeakTroughRatio = abs(pre / t);
    % trough - peak ratio

    % repol slope: the slope of the linear fit between the trough and 0.3 ms after the trough	mV/ms
    x =win_tvec(z_idx:nearest_idx(0.3,win_tvec));
    y = wave(z_idx:nearest_idx(0.3,win_tvec));
    repol_coeffs = polyfit(x, y, 1); % linear
    repol_slope = repol_coeffs(1);
    repol_int = repol_coeffs(2);

    repol_fit = polyval(repol_coeffs, x);


    % recovery slope: 	the slope of the linear fit between the post-trough peak and 0.3 ms after the peak	mV/ms
    x =win_tvec(p_idx:nearest_idx(win_tvec(p_idx)+0.3,win_tvec));
    y = wave(p_idx:nearest_idx(win_tvec(p_idx)+0.3,win_tvec));
    recov_coeffs = polyfit(x, y, 1); % linear
    recov_slope = recov_coeffs(1);
    recov_int = recov_coeffs(2);

    recov_fit = polyval(recov_coeffs, x);


    % plot
    if plot_flag

        figure(1010)
        clf
        hold on
        plot(win_tvec, wave, 'k', LineWidth=2)
        plot(win_tvec(p_idx) ,p,".c", 'MarkerSize',20)
        plot(win_tvec(t_idx) ,t,".g", 'MarkerSize',20)
        plot(win_tvec(pre_idx) ,pre,".y", 'MarkerSize',20)

        plot([win_tvec(p_idx) win_tvec(p_idx)],[p p*1.25],"-k")

        plot([win_tvec(p_idx) win_tvec(p_idx)],[p p*1.25],"-k")
        text(win_tvec(p_idx),p*1.25,"peak", 'VerticalAlignment','baseline')

        plot([win_tvec(t_idx) win_tvec(t_idx+3)],[t t],"-k")
        text(win_tvec(t_idx+3),t,"trough", 'VerticalAlignment','middle')

        plot([win_tvec(pre_idx) win_tvec(pre_idx)],[pre pre*1.25],"-k")
        text(win_tvec(pre_idx),pre*1.25,"pre-peak", 'VerticalAlignment','baseline')

        %peak trough
        plot(win_tvec([t_idx,p_idx]),[p,p],"-b")
        text(win_tvec(t_idx), p*1.1, ['dur: ' num2str(dur,2) 'ms'], 'VerticalAlignment','baseline')

        plot(win_tvec([t_idx,t_idx]),[p,t],"-b")
        text(win_tvec(t_idx), (p-abs(t))/2, ['Peak-Trough ratio: ' num2str(troughPeakRatio,2) ], 'VerticalAlignment','baseline', 'Rotation',90)

        % pre peak to trough
        plot(win_tvec([pre_idx,t_idx]),[t,t],"-m")

        plot(win_tvec([pre_idx,pre_idx]),[pre,t],"-m")
        text(win_tvec(pre_idx), (pre-abs(t))/2, ['prePeak-Trough ratio: ' num2str(prePeakTroughRatio,2) ], 'VerticalAlignment','baseline', 'Rotation',90)


        % text((win_tvec(t_idx)-win_tvec(p_idx))+win_tvec(p_idx),p*1.1,num2str(width))

        % repol
        % plot(win_tvec(z_idx:nearest_idx(0.3,win_tvec)), wave(z_idx:nearest_idx(0.3,win_tvec)), 'r')
        plot(win_tvec(z_idx:nearest_idx(0.3,win_tvec)), repol_fit, '--r', 'LineWidth',2)

        % recovery
        plot(win_tvec(p_idx:nearest_idx(win_tvec(p_idx)+0.3,win_tvec)), recov_fit, '--r', 'LineWidth',2)


        title(['Cell #: ' num2str(ii)])

            % pause
    end



    %% collect the features

    S.usr.wave(ii,:) = wave;
    S.usr.wave_tvec(ii,:) = win_tvec; 
    S.usr.peak(ii) = p;
    S.usr.trough(ii) = t;
    S.usr.pre(ii) = pre;
    S.usr.amp(ii) = amp;
    S.usr.dur(ii) = dur;
    S.usr.pre_peak_ratio(ii) = pre_peak_r; 
    S.usr.pre_t_dur(ii)  = win_tvec(t_idx) - win_tvec(pre_idx); 
    S.usr.pt_ratio(ii) = troughPeakRatio;
    S.usr.pre_t_ratio(ii) = prePeakTroughRatio;
    S.usr.repol_slope(ii) = repol_slope;
    S.usr.recov_slope(ii) = recov_slope;


    %% try some basic dimentionality reduction and clustering

    % map_dat = [dur; troughPeakRatio; prePeakTroughRatio; repol_slope; recov_slope]; 

    % wave_map = umap
end