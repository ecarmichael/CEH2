function MS_ModIdx(sig_phi, sig_amp, win_s)
%% MS_ModIdx: computes the modulation index as a continuous variable between two signals. Indended for theta - gamma modulation as per Tort et al. 2010/2009/2008
%
%
%
%    Inputs: 
%    - sig_phi: [struct] filtered LFP data in the csc format (vdmlab codebase) to be used for the phase component.  
%
%    - sig_amp: [struct] filtered LFP data in the csc format (vdmlab codebase) to be used for the amplitude component. 
%
%    - win_s: [double]   window size in samples. 
%
%    Outputs: 
%    - ModIdx: [nSamples x 1]   modulation index. 
%
%
%
%
% EC 2022-04-26   initial version 
%
%
%
%% initialize

% load some data to try out
% 
% 
% data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Inter\pv1060\HATD1'; 
% cd(data_dir)
% 
% load('LFP_mats\all_t_post_REM.mat')
% load('LFP_mats\all_MG_post_REM.mat')
% load('LFP_mats\all_t_freq_post_REM.mat')
% load('LFP_mats\all_MG_freq_post_REM.mat')

%% extract phase and amplitude

Phi_1 = angle(hilbert(sig_phi.data)); 
Amp_1 = smooth(abs(hilbert(sig_amp.data)), floor(sig_amp.cfg.hdr{1}.SamplingFrequency*0.1))'; 

%% check the signals 
plot_win = 1:(2*sig_phi.cfg.hdr{1}.SamplingFrequency); 
figure(801)
clf
subplot(3,2,1:2)
hold on
yyaxis left
plot(sig_phi.tvec(plot_win)/1000, sig_phi.data(plot_win), 'b')
plot(sig_amp.tvec(plot_win)/1000, sig_amp.data(plot_win), 'r')
plot(sig_amp.tvec(plot_win)/1000, Amp_1(plot_win), '--r')
ylabel('amplitude')

yyaxis right
plot(sig_amp.tvec(plot_win)/1000, Phi_1(plot_win), '--b')
ylim([ -2*pi 2*pi])
ylabel('phase')
xlim([sig_phi.tvec(plot_win(1)) sig_phi.tvec(plot_win(end))]/1000); 
legend({'Phi signal', 'Amp signal'})
title('filtered signals w/amplitude')

%% sliding window MI calculations
phi_bins = -pi:pi/9:pi; 

ModIdx = NaN(size(Phi_1)); 
    tic
for ii = 1:win_s:length(ModIdx)
    if (ii+win_s)> length(ModIdx)
        this_idx = ii:length(ModIdx); 
    else
        this_idx = ii:ii+win_s; 
    end
        this_phi = Phi_1(this_idx); 
        this_amp = Amp_1(this_idx); 
    
    [~,phi_idx] = histc(this_phi, phi_bins); 

    % split this block into phase bins    
    for iB = unique(phi_idx)
        phi_amp(iB) = nanmean(this_amp(phi_idx == iB)); 
    end
    
    shuff_amp = this_amp(randperm(length(this_amp))); 
    for iB = unique(phi_idx)
        phi_amp_shuff(iB) = nanmean(shuff_amp(phi_idx == iB)); 
    end
    
    % divide average amp by total over bins to normalize
%     phi_amp = phi_amp./sum(phi_amp); 
    
    nbins = length(phi_amp); 

%     Pj = (phi_amp+eps)/sum(phi_amp); % normalize phase-amp coupling
%     Qj = repmat(mean(Pj), 1, length(Pj)); % normal distribution for Pj. 
%     Hp = -sum(Pj.*log(Pj)); % shannon entropy H for P
%     Dkl = log(nbins) - Hp; % KL distance conversion of Hp
%     Dkl = sum(Pj.*log(Pj./Qj)); 
%     MI = Dkl./log(nbins);  % get modulation index. 
    
    
    ModIdx(this_idx) = (log(nbins)-(-sum((phi_amp/sum(phi_amp)).*log((phi_amp/sum(phi_amp))))))/log(nbins); 
%     ModIdx(this_idx)= (log(nbins) - (-sum(((phi_amp+eps)/sum(phi_amp)).*log((phi_amp+eps)/sum(phi_amp)))))    /log(nbins);
%     ModIdx(this_idx)= (log(nbins) - (-sum(((phi_amp_shuff+eps)/sum(phi_amp_shuff)).*log((phi_amp_shuff+eps)/sum(phi_amp_shuff)))))    /log(nbins);

%     MI_sgamma = (log(nbin)-(-sum((PAC_slow_gamma/sum(PAC_slow_gamma)).*log((PAC_slow_gamma+eps/sum(PAC_slow_gamma))))))/log(nbin); % from GE. 
end
toc
ModIdx_z = zscore(ModIdx); 
%% check the output

figure(1010)
clf
ax(1) = subplot(4,1,1:2);
hold on
plot(sig_phi.tvec, Phi_1)
plot(sig_phi.tvec, Amp_1*10000);

ax(2) = subplot(4,1,3);
plot(sig_amp.tvec, ModIdx_z)

% ax(3) = subplot(4,1,4);
% plot(sig_amp.tvec, ModIdx)

linkaxes(ax, 'x')