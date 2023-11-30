function [CoMo,Phi_f, Amp_f] = MS_phase_freq(cfg_in, csc, phi_range, amp_range)
%% MS_phase_freq:  Recreate Phase-freq Modulation IDX plot (FIG 1B) from Tort et al. 2009 https://www.pnas.org/doi/full/10.1073/pnas.0911331106
%
%   "The comodulgram plot shown in Fig. 1B was obtained by applying the MI
%   measure to several narrowed-filtered frequency pairs 
%   (phase frequencies: 4-Hz bandwidths and 2-Hz steps; amplitude 
%   frequencies: 10-Hz bandwidths and 5-Hz steps) "
%
%
%
%    Inputs: 
%    - cfg [struct]   configuration see the defaults below. 
%
%    - csc [struct]   LFP data in the TSD format
%
%    - phi_range:[min_freq max_freq]   range of phase frequencies
%
%    - amp_range:[min_freq max_freq]   range of amp frequencies
%
%    Outputs: 
%    - CoMod: [phi_range X amp_range]  MI values per phi-amp pair. 
%
%
%
%
% EC 2023-11-07   initial version 
%
%
%
%% initialize
cfg_def = [];
cfg_def.P_win = [-2 2]; % range around the Phi freq of interest. 
cfg_def.A_win = [-5 5]; % range around the Amp freq of interest. 
cfg_def.P_step = 2; % step size in Hz
cfg_def.A_step = 5; % step size in Hz
cfg_def.channel = 1; % only use the first csc channel; 
cfg_def.phi_bins = 18; 
cfg = ProcessConfig(cfg_def, cfg_in);
%%  loop over over phase ranges and amplidues 

CoMo = []; 

Phi_f = phi_range(1):cfg.P_step:phi_range(end);
Amp_f = amp_range(1):cfg.A_step:amp_range(end);


for iP = length(Phi_f):-1:1
    
    this_f = Phi_f(iP)+cfg.P_win; 
    
    if this_f(2) <=14
        method = 'cheby1';
    else
        method = 'butter';
    end
    
    
    cfg_filt = [];
    cfg_filt.type = method;%'fdesign'; %the type of filter I want to use via filterlfp
    cfg_filt.f  = this_f; % freq range to match Mizuseki et al. 2011
    cfg_filt.order = 3; %type filter order
    cfg_filt.display_filter = 0; % use this to see the fvtool
    cfg_filt.verbose = 0;
    
    csc_f = FilterLFP(cfg_filt, csc); % filter the raw LFP using
    
    
    phi = angle(hilbert(csc_f.data(cfg.channel,:)));
    
    
    % loop over amp frequencies
    
    for iA = length(Amp_f):-1:1
        
        
        
        %fprintf('Processing: Phi %.1f Amp %.1f\n',Phi_f(iP), Amp_f(iA))
        
        
        this_f = Amp_f(iA)+cfg.A_win;
        
        if this_f(2) <=14
            method = 'cheby1';
        else
            method = 'butter';
        end
        
        
        cfg_filt = [];
        cfg_filt.type = method;%'fdesign'; %the type of filter I want to use via filterlfp
        cfg_filt.f  = this_f; % freq range to match Mizuseki et al. 2011
        cfg_filt.order = 3; %type filter order
        cfg_filt.display_filter = 0; % use this to see the fvtool
        cfg_filt.verbose = 0;
        
        csc_f = FilterLFP(cfg_filt, csc); % filter the raw LFP using
        
        
        amp = abs(hilbert(csc_f.data(cfg.channel,:)));
        
        
        % compute the ModIDX
        phi_bins = -pi:pi/cfg.phi_bins:pi;
        [~,bins_idx]= histc(phi, phi_bins); %This creates a vector with the corresponding bin of each phase
        %3. The phases are binned according to which phase they belong to and their
        %mean is calculated
        amp_means= zeros(1,length(unique(bins_idx)));
        for ii= 1: length(unique(bins_idx));
            amp_means(ii)= nanmean(amp(bins_idx == ii));
        end
        %4.Each mean amplitude is then divided by the sum over the bins
        phi_amp= amp_means/sum(amp_means);
        
        
        CoMo(iP, iA) =MS_ModIdx(phi_amp);
%         amp_val(iA) = 
    end
    
    
    
    
end

%% plot if you like
hold on
imagesc(Phi_f, Amp_f, CoMo)
plot(mean(CoMo,2))
