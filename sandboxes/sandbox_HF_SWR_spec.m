%% sandbox HF_SWR_spec

csc_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\PoxR1\HF\Pox3265_2026-05-12_14-14-58_SWR2\Record Node 117'; 

swr_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\PoxR1\HF\Pox3265_2026-05-12_14-14-58_SWR2\Record Node 143'


csc_idx = [];
ts_prime = 0; 
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

%% offline SWR detection

swrs = MS_SWR_detector(csc, 'CH64')
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

vline(SWR_evts.t{2}, 'r')

%% SWR-Trigged average
csc_idx = find(contains(csc.label, 'CH57'))

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

hold
plot((0:length(STA)-1)./csc.cfg.hdr{1}.SamplingFrequency, mean(STA, 'omitmissing'))
plot((0:length(STA)-1)./csc.cfg.hdr{1}.SamplingFrequency, mean(STA, 'omitmissing') +  std(STA), '--r')
plot((0:length(STA)-1)./csc.cfg.hdr{1}.SamplingFrequency, mean(STA, 'omitmissing') -  std(STA), '--r')


xlim([.5 1.5])


%% get Single Double Tripples as per Yamamato & Tonegawa Neuron 2017
% https://www.cell.com/neuron/fulltext/S0896-6273(17)30857-7#sec-4

swr_d = IVcenters(swrs); 
single_off = swr_d > .200 | swr_d < .70; 
double_off = swr_d < .200 & swr_d > .70; 
triple_off = swr_d > .200; 


oe_swr_f = diff(SWR_evts.t{2}); 

swr_type  = NaN(size(swr_d));
for ii = 1:length(swr_d)-1
    if swr_d(ii+1) - swr_d(ii) > .200
    swr_type(ii) = 1;
    elseif ((swr_d(ii+1) - swr_d(ii)) < .200) && ((swr_d(ii+1) - swr_d(ii)) >.70)
    swr_type(ii) = 2;
    else
    swr_type(ii) = 0;
    end
end

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


