function mod_out = MS_phase_mod(cfg_in, csc, S, k_idx, plot_flag)

%% MS_phase_mod:
%
%
%
%    Inputs:
%    - cfg_in: [struct] parameters if empty uses defaults.
%    - S: [struct]  spikes in the ts format
%    - csc: [struct]  LFP data in the csc tsd format.
%
%
%
%    Outputs:
%    -

%
%
%
%
% EC 2026-08-27   initial version
%
%
%
%% initialize

if nargin < 3
    k_idx = 1:length(csc.tvec);
    plot_flag = 0;
elseif nargin < 4
    plot_flag = 0;
end

cfg_def = [];
cfg_def.f = [6 10]; % frequencies
cfg_def.type = 'butter';
cfg_def.order = 4;
cfg_def.band = 'bandpass'; % only implemented for butterworth
cfg_def.R = 0.5; % passband ripple (in dB) for Chebyshev filters only

% power threshold
cfg_def.pow_thresh = []; % if empty don't use.


cfg = ProcessConfig2(cfg_def, cfg_in);


%% Filter the LFP in the band of interest

csc_f = FilterLFP(cfg, csc);


phase = angle(hilbert(csc_f.data(1,:)));

amp = abs(hilbert(csc_f.data(1,:)));

if ~isempty(cfg.pow_thresh)
    k_idx = amp > cfg.pow_thresh;
end

%%

fprintf('\nCell #    ')

for iC = length(S.t):-1:1

    s_phi = interp1(csc_f.tvec, phase, S.t{iC});
    [s_phi_h,  b] = histcounts(s_phi, -pi:(pi/16):pi);


    mod_out.MVL(iC) = circ_r(s_phi_h,s_phi_h, -pi:(pi/16):pi); % get the mean vector length as a measure of how strong the theta locking is.

    mod_out.mean(iC) = circ_mean(s_phi);
    mod_out.phi{iC} = s_phi;
    mod_out.phi_s(iC,:) = s_phi_h;

    % grabs some shuffle values.
    for i_p_shuff = 500:-1:1
        mod_out.shuff_phi{iC}(i_p_shuff,:) = randsample(phase, length(S.t{iC}));
        mod_out.shuff_phi_h{iC}(i_p_shuff,:) = histcounts(mod_out.shuff_phi{iC}(i_p_shuff,:), -pi:(pi/16):pi);
    end

    mod_out.shuff_MVL{iC} = circ_r(mod_out.shuff_phi{iC},[],2);
    mod_out.shuff_mean{iC} = circ_mean(mod_out.shuff_phi{iC},[],2);
    mod_out.shuff_mean_h{iC} = circ_mean(mod_out.shuff_phi_h{iC},[],1);
    
    % zscore the data to the shuffle
    mod_out.MVL_z(iC) = (mod_out.MVL(iC) - mean(mod_out.shuff_MVL{iC})) / std(mod_out.shuff_MVL{iC});
    mod_out.phi_z(iC,:) = (mod_out.phi_s(iC,:) - mean(mod_out.shuff_mean_h{iC})) / std(mod_out.shuff_mean_h{iC});

    % counter
    fprintf('\b\b\b\b\b%2d/%2d', iC, length(S.t));

end
fprintf('\b\b\b\b\b%2d/%2d  - done\n', 0, length(S.t));

%%
if plot_flag

    figure(7677)
    clf
    % subplot(10,6,[1 2 7 8 13 14 19 20])
    imagesc([rad2deg(b) rad2deg(b)+180], 1:length(S.t), [normalize(mod_out.phi_z,2) normalize(mod_out.phi_z,2)])%[zscore(mod_out.phi_s, [], 2) zscore(mod_out.phi_s, [], 2)])
    % xlim([-180 180]);
    set(gca, 'XTick', [-180:90:540])
    ylabel('Cell #'); xlabel('Phase (deg)')

    figure(7676)
    clf
    m = ceil(length(S.t)/6);
    n = ceil(length(S.t)/8);

    for ii = 1:length(S.t)
        subplot(n, m, ii)
        if mod_out.MVL_z(iC) > 2
            polarplot(mod_out.phi_s(ii,:), 'LineWidth',2);
        else
            polarplot(mod_out.phi_s(ii,:));

        end
        hold on
        polarplot(mod_out.shuff_mean_h{iC});

        set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise', 'FontSize', 4)


    end


end



%% outtakes

% t_d = [];
%    for ii = length(root_baseline_stim.b_lfp):-1:1
%        d= abs(hilbert(root_baseline_stim.b_lfp(ii).DeltaFilter(root_baseline_stim.b_lfp(ii).signal, root_baseline_stim.b_lfp(ii).fs)));
%        t_d(ii) = mean(root_baseline_stim.b_lfp(ii).theta_amplitude ./d);
%    end
%    [~, td_idx] = max(t_d);
%    data.theta_phi = root_baseline_stim.b_lfp(td_idx).theta_phase;
%    data.theta_ts = root_baseline_stim.b_lfp(td_idx).ts;