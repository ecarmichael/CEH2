%% JB SWR screener


%% init

% data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\B10_chrna2_electrophy_tungtsen electrode\B10_chrna2_SMpredictible_D1\GcfB10_Chrna2_SM-PUFF_D1_915';
data_dir = '/home/williamslab/Dropbox (Williams Lab)/Data SWR sm/JB1556_2022-03-09_day1';




cd(data_dir);


%% screen the PSDs (only run this once. It will save a png of all the PSDs)

MS_Quick_psd

%% get the session information


parts = strsplit(cd, filesep);

parts = strsplit(parts{end}, '_');

f_info.subject = parts{1};
f_info.date = parts{2};
f_info.session = parts{end};



%% load the data

cfg_csc = [];
if strcmp(f_info.subject, 'JB1556')
    cfg_csc.fc = {'CSC8.ncs'};%, 'CSC3.ncs'};  % pick the channels to load
    conv_fac = 138/30;
else
    cfg_csc.fc = {'CSC16.ncs'}; % for JB 1446;
end
csc = LoadCSC(cfg_csc); % load the csc data

cfg_pos = [];
cfg_pos.convFact = [conv_fac conv_fac];
pos = LoadPos(cfg_pos); % load the position data.  This appears to be empty.

evt = LoadEvents([]);


% restrict csc and position to 4 hours.
csc = restrict(csc, csc.tvec(1), csc.tvec(1)+(4*60*60));
pos = restrict(pos, csc.tvec(1), csc.tvec(1)+(4*60*60));


% get the movement
linspeed = getLinSpd([], pos);

speed_int = interp1(linspeed.tvec, smooth(linspeed.data, csc.cfg.hdr{1}.SamplingFrequency*2), csc.tvec);

move_idx = (speed_int >1)';

%% get the theta delta ratio
cfg_con_f = [];
cfg_con_f.threshold = 0;
cfg_con_f.f = [6 12]; cfg_con_f.type = 'fdesign';
%     cfg_con_f.display_filter = 1
csc_th = FilterLFP(cfg_con_f,csc);
% csc_th.data = smooth(csc_th.data, csc.cfg.hdr{1}.SamplingFrequency*2);
csc_th.data = smooth(abs(hilbert(csc_th.data)), csc.cfg.hdr{1}.SamplingFrequency*10);

cfg_con_f = [];
cfg_con_f.threshold = 0;
cfg_con_f.f = [1 4]; cfg_con_f.type = 'fdesign';
%     cfg_con_f.display_filter = 1
csc_delta = FilterLFP(cfg_con_f,csc);
% csc_delta.data = smooth(csc_delta.data, csc.cfg.hdr{1}.SamplingFrequency*2);
csc_delta.data = smooth(abs(hilbert(csc_delta.data)), csc.cfg.hdr{1}.SamplingFrequency*10);


ratio = csc;
ratio.data = zscore((csc_th.data ./ csc_delta.data)');

z_ratio = nan(1, length(csc.tvec));
z_ratio(~move_idx) = zscore(csc_th.data(~move_idx)./csc_delta.data(~move_idx));

%% filter the LFP into the ripple band.
cfg_swr = [];
cfg_swr.check = 0; % plot checks.
cfg_swr.filt.type = 'butter'; %Cheby1 is sharper than butter
cfg_swr.filt.f  = [125 200]; % broad, could use 150-200?
cfg_swr.filt.order = 4; %type filter order (fine for this f range)
cfg_swr.filt.display_filter = 0; % use this to see the fvtool

% smoothing
cfg_swr.kernel.samples = csc.cfg.hdr{1}.SamplingFrequency/100;
cfg_swr.kernel.sd = csc.cfg.hdr{1}.SamplingFrequency/100;

% detection
% cfg_swr.artif_det.method = 'zscore';
% cfg_swr.artif_det.threshold = 8;
% cfg_swr.artif_det.dcn = '>';
% cfg_swr.artif_det.rm_len = .2;
cfg_swr.threshold = 2.5;% in sd
cfg_swr.method = 'zscore';
cfg_swr.min_len = 0.04; % mouse SWR: 40ms from Vandecasteele et al. 2014
cfg_swr.merge_thr = 0.02; %merge events that are within 20ms of each other.

sat_idx = (csc.data == max(csc.data)) | (csc.data == min(csc.data));

cfg_swr.nan_idx = [sat_idx | move_idx | z_ratio > 1]; % where are any nans, say from excluding artifacts, other events...

% restrictions
cfg_swr.max_len = [];
cfg_swr.max_len.operation = '<';
cfg_swr.max_len.threshold = .1;

%                 cfg_swr.min_len = [];
%                 cfg_swr.min_len.operation = '<';
%                 cfg_swr.min_len.threshold = .2;
cfg_swr.nCycles = 5; % number of cycles
cfg_swr.nCycles_operation = '>='; % number of cycles


% variaence
cfg_swr.var = [];
cfg_swr.var.operation = '<';
cfg_swr.var.threshold = 1;

SWR_evts = MS_get_LFP_events_sandbox(cfg_swr, csc);

%
%         cfg_max_len = [];
%         cfg_max_len.operation = '>';
%         cfg_max_len.threshold = 5;
%          SWR_evts = SelectIV(cfg_max_len,SWR_evts,'nCycles');

cfg_plot.display = 'iv';
cfg_plot.title = 'var_raw';
PlotTSDfromIV(cfg_plot, SWR_evts, csc)

keep_idx_tsd = csc;
keep_idx_tsd.data = [sat_idx | move_idx | z_ratio > 1];


%% check figure
cord = linspecer(3);
figure(102)
hold on
% yyaxis right
plot((linspeed.tvec -linspeed.tvec(1))/60/60 , linspeed.data,  'color', cord(2,:));
% ylim([0 50])
% yyaxis left
plot((csc.tvec - csc.tvec(1))/60/60, csc.data*1000, 'color', cord(1,:));
plot((csc.tvec - csc.tvec(1))/60/60, z_ratio,  'color', cord(3,:)); 
plot((csc.tvec - csc.tvec(1))/60/60, [sat_idx | move_idx | z_ratio > 1], 'k'); 
legend({ 'smooth speed','data', 'theta/delta z', 'excluded idx'}); 
xlim([min((csc.tvec - csc.tvec(1))/60/60) max((csc.tvec - csc.tvec(1))/60/60)])
xlabel('time from cp21/vehicle (hrs)')
ylabel('LFP voltage')

%% brek up into 20min blocks
nEvts = []; ndur = [];
dt = 30*60;
t = csc.tvec(1):dt:csc.tvec(end);

for ii = length(t)-1:-1:1
    if ii == length(t)-1
        these_swr = restrict(SWR_evts, t(ii), csc.tvec(end));
        these_keep = restrict(keep_idx_tsd, t(ii), csc.tvec(end));

    else
        % retrict to this time block
        these_swr = restrict(SWR_evts, t(ii), t(ii+1));
        these_keep = restrict(keep_idx_tsd, t(ii), t(ii+1));

        
        
    end
    
    ndur(ii) = sum(~these_keep.data)/csc.cfg.hdr{1}.SamplingFrequency; % convert used samples to total time used for SWR detection in this block; 
    nEvts(ii) = length(these_swr.tend); 
    
    if ndur(ii) < 120
        ndur(ii) = NaN;
        nEvts(ii) = NaN;
    end
    
%     cfg_plot = [];
%     cfg_plot.title = 'var_raw';
%     PlotTSDfromIV(cfg_plot, these_swr, csc)
    
    disp(num2str(nEvts(ii)/ndur(ii)))
    
    
    
end

figure(202)
yyaxis left
plot(t(1:end-1)/60/60, nEvts./ndur)
xlabel('time from cp21/vehicle (hrs)')
ylabel('SWR events / second')
yyaxis right
plot(t(1:end-1)/60/60, ndur)
ylabel('time used for detection per block (s)')


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

