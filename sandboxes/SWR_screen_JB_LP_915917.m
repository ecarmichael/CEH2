%% JB SWR screener


%% init

% data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\B10_chrna2_electrophy_tungtsen electrode\B10_chrna2_SMpredictible_D1\GcfB10_Chrna2_SM-PUFF_D1_915';
% data_dir = '/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/Chrna_SWR/DREADDA_chrna2E_SM_D1_1556';

% data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Data SWR sm\DREADDA_chrna2E_SM_D4(probe)_1556_electrophy';
data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Data SWR sm\DREADDA_chrna2E_SM_D3_1556_electrophy';

% data_dir = 'C:\Users\ecarm\Desktop\JB1466_2022-02-26_swr_test'
cd(data_dir);


%% screen the PSDs (only run this once. It will save a png of all the PSDs)

MS_Quick_psd


%% load the data

cfg_csc = [];
cfg_csc.fc = {'CSC10.ncs'}%, 'CSC3.ncs'};  % pick the channels to load; 
% cfg_csc.desired_sampling_frequency = 2000;
csc = MS_LoadCSC(cfg_csc); % load the csc data

% pos = LoadPos([]); % load the position data.  This appears to be empty.

evt = LoadEvents([]);

% %% high pass
% d = fdesign.highpass('N,F3dB',4,2,csc.cfg.hdr{1}.SamplingFrequency);
% Hd = design(d,'butter');
% b = Hd.sosMatrix; a = Hd.scaleValues;
% % fvtool(Hd)
%
% temp_sig = csc.data;
%
% nan_idx = find(isnan(temp_sig));
%
% if ~isempty(nan_idx)
%     fprintf('WARNING: FilterLFP.m: signal %d contains NaNs (%d).\n',iS,length(nan_idx));
%     temp_sig(nan_idx) = 0;
% end
%
% % filter
% %temp_sig = filtfilt(sos,g,temp_sig);
% temp_sig = filtfilt(b,a,temp_sig);
%
% % reinstate NaNs and put signal back into tsd
% temp_sig(nan_idx) = NaN;
% csc.data(1,:) = temp_sig;

% get the position and interp over jumps

% pos = LoadPos([]);
% 
% pos_rem = (pos.data(1,:) < 200) & (pos.data(2,:) > 375); 
% 
% plot(pos.data(1,:), pos.data(2,:), '.k')
% 
% pos_int = pos;
% pos_int.data(:,pos_rem) = NaN; 
% pos_int.data(1,:) = fillmissing(pos_int.data(1,:),'linear');
% pos_int.data(2,:) = fillmissing(pos_int.data(2,:),'linear');
% 
% hold on
% plot(pos_int.data(1,:), pos_int.data(2,:), '.r')


% linspeed = getLinSpd([],pos); % linear speed
% 
% linspeed.data = smooth(linspeed.data, ceil(1/mode(diff(linspeed.tvec)))*1)';
%  
% % Threshold speed
% cfg = []; cfg.method = 'raw'; cfg.operation = '>'; cfg.threshold = 50; % speed limit in cm/sec
% iv_fast = TSDtoIV(cfg,linspeed); % only keep intervals with speed above thresh

% interp the position, 




%% filter the LFP into the ripple band.


close all
% low pass filter the data
%     cfg_filt_d = [];
%     cfg_filt_d.type = 'butter'; %Cheby1 is sharper than butter
%     cfg_filt_d.f  = [140 250]; % broad, could use 150-200?
%     cfg_filt_d.order = 4; %type filter order (fine for this f range)
%     cfg_filt_d.display_filter = 0; % use this to see the fvtool
%     cfg.filt.display_filter = 1; % use this to see the fvtool

cfg_swr = [];
cfg_swr.check = 0; % plot checks.
cfg_swr.filt.type = 'butter'; %Cheby1 is sharper than butter
cfg_swr.filt.f  = [120 200]; % broad, could use 150-200?
cfg_swr.filt.order = 4; %type filter order (fine for this f range)
cfg_swr.filt.display_filter = 0; % use this to see the fvtool


% smoothing
cfg_swr.kernel.samples = csc.cfg.hdr{1}.SamplingFrequency/100;
cfg_swr.kernel.sd = csc.cfg.hdr{1}.SamplingFrequency/1000;

% detection
cfg_swr.artif_det.method = 'zscore';
cfg_swr.artif_det.threshold = 6;
cfg_swr.artif_det.dcn = '>';
cfg_swr.artif_det.rm_len = .1;
cfg_swr.threshold = 2;% in sd
cfg_swr.method = 'zscore';
cfg_swr.min_len = 0.04; % mouse SWR: 40ms from Vandecasteele et al. 2014
cfg_swr.merge_thr = 0.01; %merge events that are within 20ms of each other.
%

% cfg_con_f = [];
% cfg_con_f.threshold = 0;
% cfg_con_f.f = [6 12]; cfg_con_f.type = 'fdesign';
% %     cfg_con_f.display_filter = 1
% csc_con = FilterLFP(cfg_con_f,csc);
% csc_con.data = smooth(csc_con.data, csc.cfg.hdr{1}.SamplingFrequency*2); 
% csc_con.data = zscore(abs(hilbert(csc_con.data)));
% 
% cfg_con_f = [];
% cfg_con_f.threshold = 0;
% cfg_con_f.f = [1 4]; cfg_con_f.type = 'fdesign';
% %     cfg_con_f.display_filter = 1
% csc_con_d = FilterLFP(cfg_con_f,csc);
% csc_con_d.data = smooth(csc_con_d.data, csc.cfg.hdr{1}.SamplingFrequency*2); 
% csc_con_d.data = zscore(abs(hilbert(csc_con_d.data)));
% 
% 
% ratio = csc;
% ratio.data = (csc_con.data ./ csc_con_d.data)'; 
% 
% csc_con.data = smooth(csc_con.data, csc.cfg.hdr{1}.SamplingFrequency*2); 
% csc_con.data = zscore(abs(hilbert(csc_con.data)));

% %  threshold = max(csc_con.data(iChan,:))*.95;
% cfg_swr.nan_idx = linspeed.data >20;
% cfg_swr.nan_idx =  csc_con.data >= cfg_con_f.threshold;




%         sat_idx = (csc.data == max(csc.data)) | (csc.data == min(csc.data));
%         cfg_swr.nan_idx = sat_idx; % where are any nans, say from excluding artifacts, other events...
%
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

% contrast with theta
% cfg_swr.cont_filt.f = [6 10];


SWR_evts = MS_get_LFP_events_sandbox(cfg_swr, csc);


%         cfg = [];
% cfg.operation = '>';
% cfg.threshold = 5;
%         SWR_evts = SelectIV(cfg, SWR_evts, 'nCycles');





% detect theta blocks

%         cfg_t = [];
%         cfg_t.check = 0; % plot checks.
%         cfg_t.filt.type = 'Cheby1'; %Cheby1 is sharper than butter
%         cfg_t.filt.f  = [6 12]; % broad, could use 150-200?
%         cfg_t.filt.order = 4; %type filter order (fine for this f range)
%         cfg_t.filt.display_filter = 0; % use this to see the fvtool
%
%
% cfg_plot.title = 'contrast';
cfg_plot = [];
cfg_plot.display = 'iv';
cfg_plot.width = .2;
cfg_plot.title = 'var_raw';
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

