%% sandbox HF_SWR_spec
% pox1
% csc_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\PoxR1\HF\Pox3265_2026-05-12_14-14-58_SWR2\Record Node 117'; 
% swr_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\PoxR1\HF\Pox3265_2026-05-12_14-14-58_SWR2\Record Node 143'; 
% csc_idx = [24, 34];

% % pox 1_2
% phy_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_14-03-46_SWR/kilosort4'; 
% csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_14-03-46_SWR/Record Node 117'; 
% swr_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_14-03-46_SWR/Record Node 143'; 
% csc_idx = [ 15    26     28    30    32    34]; % 26

% pox_2_4
phy_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-53-28_SWR2/kilosort4'; 
csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-53-28_SWR2/Record Node 117'; 
swr_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/HF/Pox3265_2026-05-12_15-53-28_SWR2/Record Node 143'; 
swr_cur_iv = iv([csc.tvec(1) 3460 3484.5 3513.35 3805], [3400 3483.5 3512.35 3706 csc.tvec(end)]); 
theta_cur_iv = iv([3400 3706],  [3460  3805]); 
csc_idx = [ 5:16    20:30    32    33:36]; % 26

% ctrl
% csc_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\test_data\H1b32026-03-14_22-43-13_dSub_SWR\Record Node 117'

ts_prime = 0; 

%% load the spikes if present

params = OE_load_params(phy_dir);

data.S = OE_phy2TS(phy_dir, params);

ts_prime = readNPY([phy_dir filesep 'timestamps.npy']);
%% load the csc

 
if isempty(csc_dir)
    csc = []; 
    OE_evts =[]; 
else
csc_list = dir([csc_dir filesep '*CH*.continuous']);

% sort the csc based on channel number. 
for ii = length(csc_list):-1:1
    ch_idx =  strfind(csc_list(ii).name,'_CH'); 
    con_idx =  strfind(csc_list(ii).name,'.continuous'); 
    csc_num(ii) = str2double(csc_list(ii).name(ch_idx+3:con_idx)); 
    csc_name{ii} = csc_list(ii).name; 
end

% sort
[~, sort_idx] = sort(csc_num); 
csc_list = csc_list(sort_idx); 

if ~isempty(csc_idx)  && isnumeric(csc_idx(1))
    csc_list(~ismember(1:length(csc_list), csc_idx))= [];

elseif  ~isempty(csc_idx)  && iscell(csc_idx) % if csc_idx is a string look for those patterns
    
    temp_list = [];

    for ii = 1:length(csc_idx)
        temp_idx = find(contains(csc_name, [csc_idx{ii} '.continuous']));
        if length(temp_idx) > 1
            error('Too many csc files containing this name')
        end
        if ~isempty(temp_idx)
            temp_list(ii).name= csc_name{temp_idx};
            temp_list(ii).folder= csc_list(1).folder;
        end
    end

    csc_list = temp_list;

end


csc= []; labels = [];
for ii = 1:length(csc_list)

    if ii == 1
        [data, tvec, info] = load_open_ephys_data([csc_list(ii).folder filesep csc_list(ii).name]);
        csc = tsd(tvec, data);
        labels{ii} = info.header.channel;
        csc.cfg.hdr{ii} = info.header;
        csc.cfg.hdr{ii}.SamplingFrequency = info.header.sampleRate;

    else
        [data, ~, info] = load_open_ephys_data([csc_list(ii).folder filesep csc_list(ii).name]);
        csc.data =[csc.data, data];
        labels{ii} = info.header.channel;
        csc.cfg.hdr{ii} = info.header;
        csc.cfg.hdr{ii}.SamplingFrequency = info.header.sampleRate;
    end
end
fs = csc.cfg.hdr{1}.SamplingFrequency; 

csc.data = csc.data';
csc.label = labels;


% csc.tvec = ;
% csc.tvec = csc.tvec - csc.tvec(1) + (csc.tvec(1) - ts_prime(1)); % zero out the csc. 
csc.tvec = csc.tvec - ts_prime(1); % zero out the csc. 

cfg_in.decimateFactor = 15;
csc = decimate_tsd(cfg_in, csc);

% load the OE version of the events.
evts_list = dir([csc_dir filesep '*Data*.events']);

OE_evts = OE_LoadEvents([evts_list.folder filesep evts_list.name], fs);

end

% grab the SWR times
evts_list = dir([swr_dir filesep '*Data*.events']);
SWR_evts = OE_LoadEvents([evts_list.folder filesep evts_list.name], fs);

iRi = diff(SWR_evts.t{2}); 
keep_idx = iRi <.05; 

SWR_evts.t{2}(keep_idx) = []; 


if length(csc_idx) > 2
csc.data(end+1,:) = mean(csc.data); 

csc.label{end+1} = 'mean'; 
csc.cfg.hdr{end+1} = csc.cfg.hdr{end}; 
end

%% offline SWR detection

swrs = MS_SWR_detector(csc, 'CH53');

% swrs.usr.gap = [swrs.tstart(2:end) - swrs.tend(1:end-1); NaN]; 
% 
% s_idx = swrs.usr.gap
% 
% cfg_s.operation = '>=';
% cfg_s.threshold = .2; 
% swrs_s  = SelectIV(cfg_s, swrs, 'gap');
% 
% swrs_s = DifferenceIV([], swrs_s, swrs)
% 
% cfg_m.gap = .200; 
% swrs_trip = MergeIV(cfg_m, swrs); 
% 
% swrs_trip = DifferenceIV([], swrs_trip, swrs)
%% csc check

figure(1010)
clf; 
hold on
for ii = 1:size(csc.data,1)

    plot(csc.tvec, csc.data(ii,:)+ii*500);
    lab{ii} = csc.label{ii}; 
    y_t(ii) = median(csc.data(ii,:)+ii*500); 
end
set(gca, 'YTick', y_t, 'YTickLabel', lab)

% vline(SWR_evts.t{2}, 'r')

%% SWR-Trigged average
csc_idx = find(contains(csc.label, 'CH64'))

STA = []; 
win = 1; 

    startIdx = nearest_idx3(SWR_evts.t{2} - win, csc.tvec);
    endIdx = nearest_idx3(SWR_evts.t{2} + win, csc.tvec);

    

for ii = length(SWR_evts.t{2})-50:-1:50

    % Extract the segment of the csc data around the SWR event
    STA(ii, :) = csc.data(csc_idx, startIdx(ii):endIdx(ii));


end

figure(10)
clf;
subplot
hold
plot((0:length(STA)-1)./csc.cfg.hdr{1}.SamplingFrequency, mean(STA, 'omitmissing'))
plot((0:length(STA)-1)./csc.cfg.hdr{1}.SamplingFrequency, mean(STA, 'omitmissing') +  std(STA), '--r')
plot((0:length(STA)-1)./csc.cfg.hdr{1}.SamplingFrequency, mean(STA, 'omitmissing') -  std(STA), '--r')


xlim([.5 1.5])


%% get Single Double Tripples as per Yamamato & Tonegawa Neuron 2017
% https://www.cell.com/neuron/fulltext/S0896-6273(17)30857-7#sec-4
this_swr = pox2_4.swrs; 
this_csc = pox2_4.csc; 

% % swrs_c = IVcenters(this_swr); 
% % swr_type  = NaN(size(this_swr.tstart));
% % swr_d = diff(swrs_c); 
% % s_idx = swr_d > .200 & swr_d > .070; 
% % d_idx = swr_d <= .200 & swr_d > .070; 
% % 
% % swr_type(s_idx) = 1; 
% % swr_type(d_idx) = 2; 
% % 
% % t_idx = false(size(swr_type)); 
% % for ii = 1:length(swr_type)-1
% %     if (swr_type(ii)) == 2 && (swr_type(ii+1) == 2)
% %         t_idx(ii) = true;
% %     end
% % end
% % 
% % swr_type(t_idx) = 3; 
% % swr_type([find(t_idx)+1 find(t_idx)+2]) = []; 
% % 
% % d_idx = false(size(swr_type)); 
% % for ii = 1:length(swr_type)-1
% %     if (swr_type(ii)) == 2 && (swr_type(ii+1) == 1)
% %         d_idx(ii) = true;
% %     end
% % end
% % 
% % swr_type(find(d_idx)+1) = []; clf


% 
% 
swr_type  = NaN(size(this_swr.tstart));

for ii = 1:length(this_swr.tstart)-1
    d = this_swr.tstart(ii+1:end) - this_swr.tend(ii); 

    if (length(d) > 2) &&(d(1) <= .200) && (d(1) >.070) && (d(2) <= .400) && (d(2) >.070)
    swr_type(ii) = 3;
    elseif (d(1) <= .200) && (d(1) >.070)
        swr_type(ii) = 2;
    elseif d(1) > .200 && (d(1) >.070)
        swr_type(ii) = 1;
    else
    swr_type(ii) = 4;
    end
end
swr_type(end) = []; 
% 
% % remove the doubles at the end of another doubles
% t_idx = zeros(size(swr_type)); 
% for ii = 1:length(swr_type)-1
%     if swr_type(ii) == 2 && swr_type(ii+1) == 2
%         t_idx(ii) = 1; 
%     end
% end
% swr_type(logical(t_idx)) = 3; 
% swr_type([find(t_idx)+2 find(t_idx)+3]) = []; 
% 
% 
% % remove the doubles at the end of another single
% d_idx = zeros(size(swr_type)); 
% for ii = 1:length(swr_type)-1
%     if swr_type(ii) == 2 && swr_type(ii+1) == 1
%         d_idx(ii) = 1; 
%     end
% end
% 
% swr_type(logical(d_idx)) = 2; 
% swr_type(find(d_idx)+2) = []; 

% % try using ivs
% for ii = 1:length(this_swr.tstart) - 1
% 
%     if (ii < length(this_swr.tstart) - 2) && ((this_swr.tstart(ii+1) - this_swr.tend(ii)) < .200) && ((this_swr.tstart(ii+2) - this_swr.tend(ii+1)) < .200)
%         swr_c(ii:ii) = 3; 
% 
%     elseif 
%     end
% 
% end


figure(101)
clf
histogram(swr_type, 'Normalization','probability')
set(gcf,'Units','inch','position',[3 3 5 4]);
ylabel('probability'); xlabel('SWR type')

swr_types_rate = [sum(swr_type == 1) ./ (length(this_csc.tvec) / this_csc.cfg.hdr{1}.SamplingFrequency),...
    sum(swr_type == 2) ./ (length(this_csc.tvec) / this_csc.cfg.hdr{1}.SamplingFrequency), ...
        sum(swr_type == 3) ./ (length(this_csc.tvec) / this_csc.cfg.hdr{1}.SamplingFrequency),...
            sum(swr_type == 4) ./ (length(this_csc.tvec) / this_csc.cfg.hdr{1}.SamplingFrequency)];

%% histo of SWR diffs

figure(1919)
subplot(2,2,1)
histogram(diff(SWR_evts.t{2}), 0:.01:1)
vline([110, 220, 330]./1000)

subplot(2,2,3)
% pie([])

subplot(2,2,2)
histogram(diff(IVcenters(swrs)), 0:.01:1)
vline([110, 220, 330]./1000)



%% get the spectrogram


%%

this_csc = pox1_2.csc;
this_swr = pox1_2.swrs; 

csc_idx = contains(this_csc.label, 'mean');

csc_ft = this_csc; 
csc_ft.data = this_csc.data(find(csc_idx),:); 
csc_ft.label = []; 
csc_ft.label{1} = this_csc.label{find(csc_idx)}; 



swr_centers = IVcenters(this_swr); 

[csc_ft_out, TFR] = Triggered_Spec_FT(csc_ft, swr_centers, [], 80:.5:200, [], [-.5 .5]);
z_idx = nearest_idx3(0, TFR.time); 
b_idx = [nearest_idx3(-.5, TFR.time), nearest_idx3(-.2, TFR.time)];  

set(gcf,'Units','inch','position',[3 3 5 4]);
xlim([-.05 .05]); ylabel('frequency (hz)'); xlabel('time from ripple center (ms)')
set(gca,'xtick', -.05:0.025:.05); 
set(gca, 'XTickLabel', get(gca, 'XTick')*1000)

for ii = length(swr_centers):-1:1
    [~, TFR] = Triggered_Spec_FT(csc_ft, swr_centers(ii), [], 80:.5:200,  [], [-.5 .5], 0);
    psds(ii, :) = (squeeze(TFR.powspctrm(1,:,z_idx)));% - mean(squeeze(TFR.powspctrm(1,:,b_idx(1):b_idx(2))), 2, 'omitmissing')')./std(mean(squeeze(TFR.powspctrm(1,:,b_idx(1):b_idx(2))), 2, 'omitmissing')');
end


%% average PSD per ripple

% this_psd = squeeze(TFR.powspctrm(1,:,z_idx)); 

figure(1013)
clf
hold on

yMean = mean(ctrl_h2.psds, 1);
ySEM = std(ctrl_h2.psds, 0, 1) / sqrt(size(ctrl_h2.psds, 1));

shadedErrorBar(ctrl_h2.TFR.freq, yMean, ySEM, 'lineprops', '-k');


yMean = mean(pox1.psds, 1);
ySEM = std(pox1.psds, 0, 1) / sqrt(size(pox1.psds, 1));

shadedErrorBar(pox1.TFR.freq, yMean, ySEM, 'lineprops', '-r');

%% session 2
