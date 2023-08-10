function ModIdx = MS_ModIdx_win(sig_phi, sig_amp, win_s, nShuff)
%% MS_ModIdx_win: computes the modulation index as a continuous variable between two signals. Indended for theta - gamma modulation as per Tort et al. 2010/2009/2008
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
%    - nShuff: [double] 
%
%    Outputs: 
%    - ModIdx: [nSamples x 1]   modulation index. 
%
%
%
% EC 2022-04-26   initial version 
%
%
%
%% extract phase and amplitude

Phi_1 = angle(hilbert(sig_phi.data)); 
Amp_1 = smooth(abs(hilbert(sig_amp.data)), floor(sig_amp.cfg.hdr{1}.SamplingFrequency*0.1))'; 

%% check the signals 
% plot_win = 1:(2*sig_phi.cfg.hdr{1}.SamplingFrequency); 
% figure(801)
% clf
% subplot(3,2,1:2)
% hold on
% yyaxis left
% plot(sig_phi.tvec(plot_win)/1000, sig_phi.data(plot_win), 'b')
% plot(sig_amp.tvec(plot_win)/1000, sig_amp.data(plot_win), 'r')
% plot(sig_amp.tvec(plot_win)/1000, Amp_1(plot_win), '--r')
% ylabel('amplitude')
% 
% yyaxis right
% plot(sig_amp.tvec(plot_win)/1000, Phi_1(plot_win), '--b')
% ylim([ -2*pi 2*pi])
% ylabel('phase')
% xlim([sig_phi.tvec(plot_win(1)) sig_phi.tvec(plot_win(end))]/1000); 
% legend({'Phi signal', 'Amp signal'})
% title('filtered signals w/amplitude')

%% limit for testsing
% phi_bins = -pi:pi/9:pi; 
% 
%     REM_idx = [575000 597000; 9250000 9628000; 31660000 32050000; 36210000 36640000]; 
% % Phi_1  = Phi_1(REM_idx(1,1):REM_idx(2)); 
% % Amp_1  = Amp_1(REM_idx(1):REM_idx(2)); 
% 
% Phi_1 = [Phi_1(REM_idx(1,1):REM_idx(1,2)) Phi_1(REM_idx(2,1):REM_idx(2,2)) Phi_1(REM_idx(3,1):REM_idx(3,2)) Phi_1(REM_idx(4,1):REM_idx(4,2))]; 
% Amp_1 = [Amp_1(REM_idx(1,1):REM_idx(1,2)) Amp_1(REM_idx(2,1):REM_idx(2,2)) Amp_1(REM_idx(3,1):REM_idx(3,2)) Amp_1(REM_idx(4,1):REM_idx(4,2))];
% 
% 
% Ls = 1:50; 
% for iL = length(Ls):-1:1
%     disp(Ls(iL))
%     win_s = (sig_phi.cfg.hdr{1}.SamplingFrequency)*(iL);
%     
%     
%     for ii = 1:win_s:length(Phi_1)
%         if (ii+win_s)> length(Phi_1)
%             this_idx = ii:length(Phi_1);
%         else
%             this_idx = ii:ii+win_s;
%         end
%         
%             this_phi = Phi_1(this_idx);
%             this_amp = Amp_1(this_idx);
%                         
%             [~,phi_idx] = histc(this_phi, phi_bins);
%             
%             % split this block into phase bins
%             for iB = unique(phi_idx)
%                 phi_amp(iB) = nanmean(this_amp(phi_idx == iB));
%             end
%             
%             ModIdx{iL}(this_idx) = MS_ModIdx(phi_amp);
%             
%     end
%     
% end
% 
% % plot the CV over 
% for iL = 1:length(Ls)
%     CV(iL) = std(ModIdx{iL}, 1)./mean(ModIdx{iL}); 
%     MI(iL) = mean(ModIdx{iL}); 
% end
% figure(919)
% subplot(1,2,1)
% plot(Ls, MI)
% 
% subplot(1,2,2)
% plot(Ls, CV)
% xlabel('time window (s)');
% ylabel('Cumulative variance')
%% sliding window MI calculations
phi_bins = -pi:pi/9:pi; 

ModIdx = NaN(size(Phi_1)); 
    tic
%    - sig_amp: [struct] filtered LFP data in the csc format (vdmlab codebase) to be used for the amplitude component. 
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
    
    
    

%     nbins = length(phi_amp); 

%     Pj = (phi_amp+eps)/sum(phi_amp); % normalize phase-amp coupling
%     Qj = repmat(mean(Pj), 1, length(Pj)); % normal distribution for Pj. 
%     Hp = -sum(Pj.*log(Pj)); % shannon entropy H for P
%     Dkl = log(nbins) - Hp; % KL distance conversion of Hp
%     Dkl = sum(Pj.*log(Pj./Qj)); 
%     MI = Dkl./log(nbins);  % get modulation index. 
    
    
    ModIdx(this_idx) = MS_ModIdx(phi_amp); % condensed formula for Mod Index here. 
    

end
toc


% %% get the shuffle values
% 
%     for iShuff  = nShuff:-1:1; 
%     shuff_amp = this_amp(randperm(length(this_amp))); 
%     for iB = unique(phi_idx)
%         phi_amp_shuff(iB) = nanmean(shuff_amp(phi_idx == iB)); 
%     end
%     end
    
%% check the output
% 
figure
clf
ax(1) = subplot(3,1,1:2);
yyaxis left
plot(sig_phi.tvec, Phi_1)
ylabel('radians')
 yyaxis right
plot(sig_phi.tvec, Amp_1*10000);
ylabel('amplitude')
title('Phase and amplitude')
legend({'Phi', 'Amp'})
set(gca, 'XTick', [])
% 
% ax(2) = subplot(4,1,3);
% plot(sig_amp.tvec, ModIdx_z)

ax(3) = subplot(3,1,3);
plot(sig_amp.tvec, ModIdx)
ylabel('Mod Idx')
legend('Mod Idx')
% 
linkaxes(ax, 'x')

xlim([sig_phi.tvec(1) sig_phi.tvec(end)])