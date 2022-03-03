%% JB SWR screener


%% init

% data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\B10_chrna2_electrophy_tungtsen electrode\B10_chrna2_SMpredictible_D1\GcfB10_Chrna2_SM-PUFF_D1_915';
data_dir = '/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/Chrna_SWR/DREADDA_chrna2E_SM_D1_1556'

data_dir = 'C:\Users\ecarm\Desktop\JB1466_2022-02-26_swr_test'
cd(data_dir);


%% screen the PSDs (only run this once. It will save a png of all the PSDs)

MS_Quick_psd


%% load the data

cfg_csc = [];
cfg_csc.fc = {'CSC15.ncs'}%, 'CSC3.ncs'};  % pick the channels to load
csc = LoadCSC(cfg_csc); % load the csc data

% pos = LoadPos([]); % load the position data.  This appears to be empty. 

evt = LoadEvents([]);

%% high pass
d = fdesign.highpass('N,F3dB',4,2,csc.cfg.hdr{1}.SamplingFrequency);
Hd = design(d,'butter');
b = Hd.sosMatrix; a = Hd.scaleValues;
% fvtool(Hd)

temp_sig = csc.data;

nan_idx = find(isnan(temp_sig));

if ~isempty(nan_idx)
    fprintf('WARNING: FilterLFP.m: signal %d contains NaNs (%d).\n',iS,length(nan_idx));
    temp_sig(nan_idx) = 0;
end

% filter
%temp_sig = filtfilt(sos,g,temp_sig);
temp_sig = filtfilt(b,a,temp_sig);

% reinstate NaNs and put signal back into tsd
temp_sig(nan_idx) = NaN;
csc.data(1,:) = temp_sig;


%% filter the LFP into the ripple band. 



% low pass filter the data
    cfg_filt_d = [];
    cfg_filt_d.type = 'butter'; %Cheby1 is sharper than butter
    cfg_filt_d.f  = [140 250]; % broad, could use 150-200?
    cfg_filt_d.order = 4; %type filter order (fine for this f range)
    cfg_filt_d.display_filter = 0; % use this to see the fvtool
    cfg.filt.display_filter = 1; % use this to see the fvtool




        cfg_swr = [];
        cfg_swr.check = 0; % plot checks.
        cfg_swr.filt.type = 'butter'; %Cheby1 is sharper than butter
        cfg_swr.filt.f  = [125 250]; % broad, could use 150-200?
        cfg_swr.filt.order = 4; %type filter order (fine for this f range)
        cfg_swr.filt.display_filter = 0; % use this to see the fvtool
        
        
        % smoothing
        cfg_swr.kernel.samples = csc.cfg.hdr{1}.SamplingFrequency/100;
        cfg_swr.kernel.sd = csc.cfg.hdr{1}.SamplingFrequency/100;
        
        % detection
        cfg_swr.artif_det.method = 'zscore';
        cfg_swr.artif_det.threshold = 4; 
        cfg_swr.artif_det.dcn = '>';
        cfg_swr.artif_det.rm_len = .2;
        cfg_swr.threshold = 2.5;% in sd
        cfg_swr.method = 'zscore';
        cfg_swr.min_len = 0.04; % mouse SWR: 40ms from Vandecasteele et al. 2014
        cfg_swr.merge_thr = 0.01; %merge events that are within 20ms of each other.
        
        sat_idx = (csc.data == max(csc.data)) | (csc.data == min(csc.data));
        cfg_swr.nan_idx = sat_idx; % where are any nans, say from excluding artifacts, other events...
        
        % restrictions
        cfg_swr.max_len = [];
        cfg_swr.max_len.operation = '<';
        cfg_swr.max_len.threshold = .2;
        
        %                 cfg_swr.min_len = [];
        %                 cfg_swr.min_len.operation = '<';
        %                 cfg_swr.min_len.threshold = .2;
        cfg_swr.nCycles = 20; % number of cycles
        cfg_swr.nCycles_operation = '<='; % number of cycles
        
        % variaence
        cfg_swr.var = [];
        cfg_swr.var.operation = '<';
        cfg_swr.var.threshold = 1;
        
        SWR_evts = MS_get_LFP_events_sandbox(cfg_swr, csc);
        
        cfg_plot.display = 'iv';
        PlotTSDfromIV(cfg_plot, SWR_evts, csc)

%% visualize
c_ord = linspecer(size(csc.data, 1)+2); 
figure(101)
hold on
for ii = 1:size(csc.data, 1)
    ax(ii) = subplot(size(csc.data, 1), 1,ii);
    hold on
plot(csc.tvec, csc.data(ii,:), 'color', c_ord(ii,:));
xlim([csc.tvec(1) csc.tvec(end)])
ylim([min(csc.data,[],'all')*1.10, max(csc.data, [], 'all')*1.10]);
plot(evt.t{3}, zeros(size(evt.t{3})), 'xk');

end

linkaxes(ax, 'x')

