%%

function OLM_SWR_screener(block_dur)
%% init

% data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\B10_chrna2_electrophy_tungtsen electrode\B10_chrna2_SMpredictible_D1\GcfB10_Chrna2_SM-PUFF_D1_915';
% data_dir = '/home/williamslab/Dropbox (Williams Lab)/Data SWR sm/JB1556_2022-03-09_day1';

% data_dir = '/home/williamslab/Dropbox (Williams Lab)/Data SWR sm/JB1556_2022-03-10_day2';

% cd(data_dir);

cord = linspecer(3);

parts = strsplit(cd, filesep);
parent_path = strjoin(parts(1:end-1), filesep);

inter_dir = [parent_path filesep 'inter'];
mkdir(inter_dir)
%% screen the PSDs (only run this once. It will save a png of all the PSDs)

% MS_Quick_psd

%% get the session information


parts = strsplit(cd, filesep);
partent = parts{1:end-1};

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
    cfg_csc.fc = {'CSC13.ncs'}; % for JB 1446;
    conv_fac = 255/30;
end
csc = LoadCSC(cfg_csc); % load the csc data

cfg_pos = [];
cfg_pos.convFact = [conv_fac conv_fac];
pos = LoadPos(cfg_pos); % load the position data.  This appears to be empty.

evt = LoadEvents([]);
start_rec = find(contains(evt.label, 'Starting Recording'));
stop_rec = find(contains(evt.label, 'Stopping Recording'));

% get the movement
linspeed = getLinSpd([], pos);

speed_int = interp1(linspeed.tvec, smooth(linspeed.data, csc.cfg.hdr{1}.SamplingFrequency*2), csc.tvec);


move_idx = (speed_int >1)';

if strcmp(f_info.subject, 'JB1556') && strcmp(f_info.session, 'day5')
    
    move_idx(nearest_idx3(18210, csc.tvec):nearest_idx3(19140, csc.tvec)) = 0; 
    move_idx(nearest_idx3(20120, csc.tvec):end) = 0;

end

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


sat_idx = (csc.data == max(csc.data)) | (csc.data == min(csc.data));


keep_idx_tsd = csc;
keep_idx_tsd.data = [sat_idx | move_idx | z_ratio > 1];
%% check figure
figure(102)
tic
hold on
% yyaxis right
plot((csc.tvec - csc.tvec(1))/60/60,  speed_int,  'color', cord(2,:));
% ylim([0 50])
% yyaxis left
plot((csc.tvec - csc.tvec(1))/60/60, csc.data*1000, 'color', cord(1,:));
plot((csc.tvec - csc.tvec(1))/60/60, z_ratio,  'color', cord(3,:));
plot((csc.tvec - csc.tvec(1))/60/60, [sat_idx | move_idx | z_ratio > 1], 'k');
legend({ 'smooth speed','data', 'theta/delta z', 'excluded idx'});
xlim([min((csc.tvec - csc.tvec(1))/60/60) max((csc.tvec - csc.tvec(1))/60/60)])
xlabel('time from cp21/vehicle (hrs)')
ylabel('LFP voltage')
toc
SetFigure([], gcf);
% saveas(gcf,[inter_dir filesep 'data_tvec'], 'png')
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
% PlotTSDfromIV(cfg_plot, SWR_evts, csc)


close all
%% restrict if needed

if length(evt.t{start_rec}) > 1
    pre_t = [evt.t{start_rec}(1), evt.t{stop_rec}(1)];
    post_t = [evt.t{start_rec}(2), evt.t{stop_rec}(2)];
    
    
    csc_pre = restrict(csc,pre_t(1), pre_t(1)+(1*60 *60));
    csc_post = restrict(csc,post_t(1), post_t(1)+(4*60 *60));
    fprintf('Restriction check: ppre duration = %0.2f hours\n', (csc_pre.tvec(end)-csc_pre.tvec(1))/60/60)
    
else
    % restrict csc and position to 4 hours.
    csc_post = restrict(csc, csc.tvec(1), csc.tvec(1)+(4*60*60));
end


fprintf('Restriction check: post duration = %0.2f hours\n', (csc_post.tvec(end)-csc_post.tvec(1))/60/60)
%% brek up into 20min blocks
nEvts = []; ndur = [];

dt = block_dur*60; % block duration in seconds.

t = csc_post.tvec(1):dt:csc_post.tvec(end)+dt; % get the time blocks in the recording.
t_zero = ((t - t(1))/60/60);%+(.5*(dt/60/60)); % zero out for plotting

t_zero = t_zero(1:nearest_idx3(4, t_zero));
t = t(1:nearest_idx3(4, t_zero));

for ii = length(t):-1:1
    if ii == length(t)
        these_swr = restrict(SWR_evts, t(ii), csc_post.tvec(end));
        these_keep = restrict(keep_idx_tsd, t(ii), csc_post.tvec(end));
    else
        % retrict to this time block
        these_swr = restrict(SWR_evts, t(ii), t(ii+1));
        these_keep = restrict(keep_idx_tsd, t(ii), t(ii+1));
    end
    
    ndur(ii) = sum(~these_keep.data)/csc.cfg.hdr{1}.SamplingFrequency; % convert used samples to total time used for SWR detection in this block;
    nEvts(ii) = length(these_swr.tend);
    
    if ndur(ii) < 300
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
plot(t_zero, nEvts./ndur)
xlabel('time from cp21/vehicle (hrs)')
ylabel('SWR events / second')
yyaxis right
plot(t_zero, ndur)
ylim([0 max(ndur)])
ylabel('time used for detection per block (s)')
xlim([0 4])
title([f_info.subject ' ' f_info.session])

saveas(gcf, [inter_dir filesep f_info.subject '_' f_info.session '_post'], 'fig')
saveas(gcf, [inter_dir filesep f_info.subject '_' f_info.session '_post'], 'png')
figure(200)
% yyaxis left
hold on
bar(t_zero-(.5*(dt/60/60)), nEvts./ndur, 'facecolor', cord(1,:), 'EdgeColor', cord(1,:))
xlabel('time from cp21/vehicle (hrs)')
ylabel('SWR rate ')

xlim([-.5 4])
plot(t_zero(isnan(nEvts)) -(.5*(dt/60/60)), zeros(1, length(t_zero(isnan(nEvts))))+max(nEvts./ndur), 'x')
title([f_info.subject ' ' f_info.session])
saveas(gcf, [inter_dir filesep f_info.subject '_' f_info.session '_post_bar'], 'fig')
saveas(gcf, [inter_dir filesep f_info.subject '_' f_info.session '_post_bar'], 'png')


% if there is a pre recording get the measures and append plot
if length(evt.t{start_rec})>1 % if there is a pre recording get it here.
    t_minus = csc_pre.tvec(1):dt:csc_pre.tvec(end);
    
    pre_nEvts = []; pre_ndur = [];
    
    t_minus_zero = ((t_minus - t_minus(1))/60/60)-1;
    
    
    
    for ii = length(t_minus):-1:1
        if ii == length(t_minus)
            these_swr = restrict(SWR_evts, t_minus(ii), csc_post.tvec(end));
            these_keep = restrict(keep_idx_tsd, t_minus(ii), csc_post.tvec(end));
        else
            % retrict to this time block
            these_swr = restrict(SWR_evts, t_minus(ii), t_minus(ii+1));
            these_keep = restrict(keep_idx_tsd, t_minus(ii), t_minus(ii+1));
        end
        
        pre_ndur(ii) = sum(~these_keep.data)/csc.cfg.hdr{1}.SamplingFrequency; % convert used samples to total time used for SWR detection in this block;
        pre_nEvts(ii) = length(these_swr.tend);
        
        if pre_ndur(ii) < 200
            pre_ndur(ii) = NaN;
            pre_nEvts(ii) = NaN;
        end
        
        disp(num2str(pre_nEvts(ii)/pre_ndur(ii)))
    end
    
    
    figure(202)
    clf
    yyaxis left
    plot([t_minus_zero    t_zero],[pre_nEvts./pre_ndur,  nEvts./ndur])
    xlabel('time from cp21/vehicle (hrs)')
    ylabel('SWR events / second')
    yyaxis right
    plot([t_minus_zero    t_zero],[pre_ndur,  ndur])
    ylim([0 max([pre_ndur, ndur])+1000])
    ylabel('time used for detection per block (s)')
    xlim([-1 4])
    rectangle('Position', [t_minus_zero(end), 0, abs(t_zero(1) - t_minus_zero(end)), max([pre_ndur, ndur])], 'facecolor', [.2 .2 .2 .2], 'edgecolor', [.2 .2 .2 .2])
    title([f_info.subject ' ' f_info.session])
    
    figure(203)
    % yyaxis left
    hold on
    bar(t_zero- (.5*(dt/60/60)), (nEvts./ndur)./nanmean(pre_nEvts./pre_ndur), 'facecolor', cord(1,:), 'EdgeColor', cord(1,:))
    xlabel('time from cp21/vehicle (hrs)')
    ylabel('SWR rate normalized to 1hr pre')
    % yyaxis right
    % plot(t_zero,  ndur)
    ylabel('time used for detection per block (s)')
    xlim([0 4])
    plot(t_zero(isnan(nEvts)) -(.5*(dt/60/60)), zeros(1,length(t_zero(isnan(nEvts))))+max(nEvts./ndur), 'x')
    title([f_info.subject ' ' f_info.session])
    
    % rectangle('Position', [t_minus_zero(end), 0, abs(t_zero(1) - t_minus_zero(end)), max([pre_ndur, ndur])], 'facecolor', [.2 .2 .2 .2], 'edgecolor', [.2 .2 .2 .2])
    saveas(gcf, [inter_dir filesep f_info.subject '_' f_info.session '_full_bar'], 'fig')
    saveas(gcf, [inter_dir filesep f_info.subject '_' f_info.session '_full_bar'], 'png')
    
end




%% collect the data for comparisons


data_out = [];
data_out.f_info = f_info;
data_out.SWR_evts = SWR_evts;
data_out.csc = csc;
data_out.keep_idx_tsd = keep_idx_tsd; 
data_out.evt = evt;
data_out.t = t;
data_out.t_zero = t_zero;



data_out.nEvts = nEvts;
data_out.nDur = ndur;

if exist('t_minus')
    data_out.pre = pre_t;
data_out.post = post_t;
    data_out.t_minus = t_minus;
    data_out.pre_nEvts = pre_nEvts;
    data_out.pre_nDur = pre_ndur;
end

save([inter_dir filesep f_info.subject '_' f_info.session '_data_out.mat'], '-v7.3', 'data_out')

close all
end