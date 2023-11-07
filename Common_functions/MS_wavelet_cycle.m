function wave_mat = MS_wavelet_cycle(cfg_in, csc); 
%% MS_wavelet_cycle:
%
%
%
%    Inputs: 
%    - cfg [struct]   configuration see the defaults below. 
%
%    - csc: [struct]: contains LFp and time data in the TSD format. 
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2023-11-03   initial version 
%
%
%
%% initialize
cfg_def = [];
cfg_def.f = [5 12]; 
cfg_def.phi_bins = 16;
cfg_def.channel = 1; 

cfg = ProcessConfig(cfg_def, cfg_in);

if cfg.f(2) <14
   cfg.method = 'cheby1';
else
   cfg.method = 'butter'; 
end
%% get wavelet for all data



[cwave, F,~,~,~] = cwt(csc.data(cfg.channel,:), csc.cfg.hdr{1}.SamplingFrequency); 

% abs_cwave = abs(cwave); % convert to magnitude. 

   
%                 AX = gca;
        [minf,maxf] = cwtfreqbounds(numel(csc.data), csc.cfg.hdr{1}.SamplingFrequency);
        
        freq = 2.^(round(log2(minf)):round(log2(maxf)));
%         AX.YTickLabelMode = 'auto';
%         AX.YTick = freq;

%% get the phase bins and create a 
cfg_filt_t = [];
cfg_filt_t.type = cfg.method;%'fdesign'; %the type of filter I want to use via filterlfp
cfg_filt_t.f  = cfg.f; % freq range to match Mizuseki et al. 2011
cfg_filt_t.order = 3; %type filter order
cfg_filt_t.display_filter = 0; % use this to see the fvtool
csc_f = FilterLFP(cfg_filt_t, csc); % filter the raw LFP using


phi = angle(hilbert(csc_f.data(cfg.channel,:)));

%% get the mean values across phase bins

phi_bins = -pi:pi/cfg.phi_bins:pi; 
[~,bins_idx]= histc(phi, phi_bins); %This creates a vector with the corresponding bin of each phase


wave_mat= nan(size(cwave,1),length(unique(bins_idx)));

for ii= 1:length(unique(bins_idx))
    wave_mat(:,ii)= nanmean(abs(cwave(:,bins_idx == ii)),2);
end


%% plot
figure(101)
clf

imagesc([phi_bins phi_bins], F, [wave_mat wave_mat])
xlabel('Phase')
ylabel('Frequency (Hz)')

AX = gca;
AX.YScale = 'log';
AX.YDir = 'normal';

caxis([ 0 0.00005])

% ylim([1 100])
% AX.YTickLabel = AX.YTick/10;

