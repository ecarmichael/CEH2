% HF_preprocessing

% check for outstanding data to preprocess

fnames = dir('/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Pox/*'); 
fnames(1:3) = []; 
for ii = length(fnames):-1:1
    Name{ii} = fnames(ii).name; 
    Path{ii} = [fnames(ii).folder filesep fnames(ii).name]; 

    if exist([Path{ii} filesep 'kilosort4'], 'dir')
        Cut(ii) = true;
    else
        Cut(ii) = false; 
    end

    if exist([Path{ii} filesep 'kilosort4' filesep 'timestamps.npy'], 'file')
        Timestamps(ii) = true;
    else
        Timestamps(ii) = false; 
    end

    % session info
    sess = Name{ii};
    if contains(sess, 'Pox')
        subject{ii} = lower(sess(strfind(sess, 'Pox'): strfind(sess, 'Pox')+6));
    elseif contains(lower(sess), 'som')
        subject{ii} = lower(sess(strfind(lower(sess), 'som'): strfind(lower(sess), 'som')+3));
    end
    if contains(sess, 'TFC')
        sess_id{ii} = strrep(sess(strfind(sess, 'TFC'):end), '_', '');
    else
        sess_id{ii} = sess(strfind(sess, 'LT'):end);
    end

    % check if the data has already been processed
    pnames = dir('/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Inter_data/HF_Sub_SWR/inter_offset/*.mat'); 

    pnames = {pnames.name}; 

    if sum(contains(lower(pnames), lower(subject{ii})) & contains(lower(pnames), lower(sess_id{ii}))) >0
        Processed(ii) = true;
    else
        Processed(ii) = false;
    end

fprintf('%s :  Kilo: %d   |   TS: %d   | Processed:  %d\n', Name{ii}, Cut(ii), Timestamps(ii), Processed(ii))
    
end

tbl = table(Name, Path, Cut, Timestamps, Processed,'VariableNames',{'Name', 'Dir', 'Cut', 'Timestamps', 'Processed'}); 


%% setup the data
if ismac
    data_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Pox/';
    inter_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Inter_data/HF_Sub_SWR/inter_offset';
elseif ispc
    data_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\Pox\';
    inter_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\Inter_data\HF_Sub_SWR\inter_offset';

end
node = [filesep 'Record Node 117'];

all_csc_dir = [];  ca1_swr_idx = []; sub_swr_idx = []; 
all_csc_dir{1} = 'Pox2217_2026-06-16_16-29-44_TFC_D1'; ca1_swr_idx{1} = {'CH59'}; sub_swr_idx{1} = {[]};
all_csc_dir{2} = 'Pox2217_2026-06-17_13-41-15_TFC_D2'; ca1_swr_idx{2} = {'CH13'}; sub_swr_idx{2} = {[]};
all_csc_dir{3} = 'Pox2217_2026-06-18_14-26-43_TFC_D3'; ca1_swr_idx{3} = {'CH119'}; sub_swr_idx{3} = {'CH157'};
all_csc_dir{4} = 'Pox2217_2026-06-19_13-37-05_TFCD4'; ca1_swr_idx{4} = {'CH59'};  sub_swr_idx{4} = {'CH157'};
all_csc_dir{5} = 'Pox2217_2026-06-20_22-52-37_TFCD5'; ca1_swr_idx{5} = {'CH124'};  sub_swr_idx{5} = {'CH138'};
all_csc_dir{6} = 'Pox3265_2026-06-16_18-09-08_TFCD1'; ca1_swr_idx{6} = {'CH51'};  sub_swr_idx{6} = {'CH141'};
all_csc_dir{7} = 'Pox3265_2026-06-17_15-59-53_TFCD2'; ca1_swr_idx{7} = {'CH63'};  sub_swr_idx{7} = {'CH133'};
all_csc_dir{8} = 'Pox3265_2026-06-18_16-10-31_TFCD3'; ca1_swr_idx{8} = {'CH51'};  sub_swr_idx{8} = {'CH141'};
all_csc_dir{9} = 'Pox3265_2026-06-19_21-06-37_TFCD4'; ca1_swr_idx{9} = {'CH63'};  sub_swr_idx{9} = {'CH143' 'CH130' 'CH142'};
all_csc_dir{10} = 'Pox3265_2026-06-20_21-00-10_TFCD5'; ca1_swr_idx{10} = {'CH56'};  sub_swr_idx{10} = {'CH143' 'CH144'};
all_csc_dir{11} = 'Pox3567_2026-06-20_17-44-01_TFCD1'; ca1_swr_idx{11} = {'CH65'};  sub_swr_idx{11} = {'CH136' 'CH134' 'CH144' 'CH141'};
all_csc_dir{12} = 'Pox3567_2026-06-21_15-41-32_TFCD2'; ca1_swr_idx{12} = {'CH124'};  sub_swr_idx{12} = {'CH149'};
all_csc_dir{13} = 'Pox3567_2026-06-22_13-44-28_TFCD3'; ca1_swr_idx{13} = {'CH74'};  sub_swr_idx{13} = {[]}; 
all_csc_dir{14} = 'Pox3567_2026-06-23_15-53-07_TFCD4'; ca1_swr_idx{14} = {'CH51'};  sub_swr_idx{14} =  {'CH142'};
all_csc_dir{15} = 'Pox3567_2026-06-24_14-10-03_TFCD5'; ca1_swr_idx{15} = {'CH71'};  sub_swr_idx{15} = {'CH134'};
all_csc_dir{16} = 'Pox3568_2026-06-20_12-49-38_TFCD1'; ca1_swr_idx{16} = {'CH51'};  sub_swr_idx{16} = {'CH143'};
all_csc_dir{17} = 'Pox3568_2026-06-21_17-05-21_TFCD2'; ca1_swr_idx{17} = {'CH119'};  sub_swr_idx{17} = {'CH157'};
all_csc_dir{18} = 'Pox3568_2026-06-22_14-58-53_TFCD3'; ca1_swr_idx{18} = {'CH118'};  sub_swr_idx{18} = {'CH145'};
all_csc_dir{19} = 'Pox3568_2026-06-23_18-12-36_TFCD4'; ca1_swr_idx{19} = {'CH70'};  sub_swr_idx{19} = {'CH145'};
all_csc_dir{20} = 'Pox3568_2026-06-24_16-56-18_TFCD5'; ca1_swr_idx{20} = {'CH75'  'CH71' 'CH115'};  sub_swr_idx{20} = {'CH149'}; % CH71 and 115 are also decent for CA1
all_csc_dir{21} = 'SOM1_2026-10-05_15-20-59_LT1';      ca1_swr_idx{21} = {''};  sub_swr_idx{21} = {[]};
all_csc_dir{22} = 'SOM1_2026-10-06_13-42-35_LT2';      ca1_swr_idx{22} = {'CH77'};  sub_swr_idx{22} = {[]};
all_csc_dir{23} = 'SOM1_2026-10-07_16-35-07_LT3';      ca1_swr_idx{23} = {''};  sub_swr_idx{23} = {[]};
all_csc_dir{end+1} = 'SOM1_2026-10-08_16-02-31_LT4';       ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'SOM1_2026-10-09_15-02-08_LT5';       ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'Som1_2026-10-12_14-35-53_TFCD1';     ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'SOM1_2026-10-13_09-36-51_TFCD2';     ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'SOM1_2026-10-14_14-47-04_TFCD3';     ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'SOM1_2026-10-15_11-42-39_TFCD4';     ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'Som1_2026-10-16_15-32-52_TFCD5';     ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
% SOM 2
all_csc_dir{end+1} = 'SOM2_2026-10-06_11-23-25_LT1';       ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'SOM2_2026-10-07_10-33-39_LT2';       ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'SOM2_2026-10-08_14-38-55_LT3';       ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'SOM2_2026-10-09_13-23-44_LT4';       ca1_swr_idx{end+1} = {'CH72'};  sub_swr_idx{end+1} = {'CH152'};
all_csc_dir{end+1} = 'SOM2_2026-10-10_12-08-36_LT5';       ca1_swr_idx{end+1} = {'CH88' 'CH102', 'CH122', 'CH116'};  sub_swr_idx{end+1} = {'CH139', 'CH132', 'CH149', 'CH155'};
all_csc_dir{end+1} = 'Som2_2026-10-12_17-07-42_TFCD1';     ca1_swr_idx{end+1} = {'CH70', 'CH75', 'CH80'};  sub_swr_idx{end+1} = {'CH156', 'CH158', 'CH160'};
all_csc_dir{end+1} = 'Som2_2026-10-13_11-41-22_TFCD2';     ca1_swr_idx{end+1} = {'CH65' 'CH77', 'CH70', 'CH76'};  sub_swr_idx{end+1} = {'CH151', 'CH145', 'CH158'};
all_csc_dir{end+1} = 'SOM2_2026-10-14_16-19-06_TFCD3';     ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'SOM2_2026-10-15_14-16-32_TFCD4';     ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'Som2_2026-10-16_19-02-55_TFCD5';     ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
% Pox3456
all_csc_dir{end+1} = 'Pox3564_2026-10-12_11-29-35_LT1';    ca1_swr_idx{end+1} = {'CH66' 'CH74', 'CH70', 'CH71', 'CH77'};  sub_swr_idx{end+1} = {'CH157', 'CH147', 'CH153', 'CH150'};
all_csc_dir{end+1} = 'Pox3564_2026-10-13_14-02-56_LT2';    ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'Pox3564_2026-10-14_09-47-57_LT3';    ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'Pox3564_2026-10-15_16-29-59_LT4';    ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'Pox3564_2026-10-16_08-43-36_LT5';    ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'Pox3564_2026-10-19_11-00-11_TFCD1';  ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
% Pox3576
all_csc_dir{end+1} = 'Pox3576_2026-10-12_12-53-56_LT1';    ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'Pox3576_2026-10-13_15-34-57_LT2';    ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'Pox3576_2026-10-14_11-26-23_LT3';    ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'Pox3576_2026-10-15_18-00-51_LT4';    ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};
all_csc_dir{end+1} = 'Pox3576_2026-10-16_13-14-05_LT5';    ca1_swr_idx{end+1} = {''};  sub_swr_idx{end+1} = {[]};


% all_csc_dir{23} = 'SOM2_2026-10-09_13-23-44_LT4';      ca1_swr_idx{23} = {'CH72'};  sub_swr_idx{23} = {'CH152'};
% all_csc_dir{24} = 'SOM2_2026-10-10_12-08-36_LT5';      ca1_swr_idx{24} = {'CH88' 'CH102', 'CH122', 'CH116'};  sub_swr_idx{24} = {'CH139', 'CH132', 'CH149', 'CH155'};
% all_csc_dir{25} = 'Pox3564_2026-10-12_11-29-35_LT1';   ca1_swr_idx{25} = {'CH66' 'CH74', 'CH70', 'CH71', 'CH77'};  sub_swr_idx{24} = {'CH157', 'CH147', 'CH153', 'CH150'};
% all_csc_dir{26} = 'Som2_2026-10-13_11-41-22_TFCD2';    ca1_swr_idx{26} = {'CH65' 'CH77', 'CH70', 'CH76'};  sub_swr_idx{24} = {'CH151', 'CH145', 'CH158'};
 

csc_idx = 1:96;

probe_type{1} = {'pox2217' 'pox3265' 'pox3567' 'pox3568'}; % A4x16 & Buz32
probe_type{2} = {'som1' 'som2' 'pox3564', 'pox3567B'};  % 5x12 & Buz32

to_cut = []; 

%% loop


for iS = 40 %[31 32 33 38 39 40]

    csc_dir = [data_dir all_csc_dir{iS} filesep 'Record Node 117'];
    swr_idx = ca1_swr_idx{iS};


    %% spikes

    if ~exist('phy_dir', 'var') && exist('csc_dir', 'var')
        if ~contains(csc_dir, '/')
            parts = strsplit(csc_dir, '\');
        else
            parts = strsplit(csc_dir, '/');
        end

        if contains(parts{1}, 'C:')
            parts(1) = [];
        end
        phy_dir = fullfile(parts{1:end-1}, 'kilosort4'); % Construct phy_dir from csc_dir
        if ispc;  phy_dir = ['C:\' phy_dir];
        elseif ismac;   phy_dir = ['/' phy_dir] ;end
    end

    % check if the session has been cut
    if  ~exist([phy_dir filesep 'cluster_info.tsv'], 'file')
        to_cut{end+1} = all_csc_dir{iS}; 
        clearvars -except all_csc_dir iS data_dir inter_dir node probe_type csc_idx ca1_swr_idx to_cut

        continue
    end

    cd(phy_dir)
    params = OE_load_params(phy_dir);

    S = OE_phy2TS(phy_dir, params);

    ch = [];
    ts_prime = readNPY([phy_dir filesep 'timestamps.npy']);
    ts_prime = ts_prime(1);

    S.loc = logical(S.usr.ch < 64);
    % c_ord = MS_linspecer(2);
    % S.c_ord = [winter(sum(S.loc==1)); flipud(nebula(sum(S.loc==0)))]%repmat(c_ord(1,:),sum(S.loc==0),1)];

    neb = nebula(round(sum(S.loc==0)*1.5));
    win = winter(round(sum(S.loc==1)*1.5));
    S.c_ord = [win(1:sum(S.loc==1),:); flipud(neb(end-sum(S.loc==0)+1:end,:))];
    %

    cd(inter_dir)

    %% load the csc

    % sort the channels based on the probe config
    % csc_names = {'CH54', 'CH4', 'CH60', 'CH58', 'CH62', 'CH15', 'CH64', 'CH56', 'CH6', 'CH11', 'CH8', 'CH52', 'CH16', 'CH12', 'CH50',...

    T = readtable([csc_dir filesep 'settings.xml']);
    csc_names = T.index_1Attribute(3,:);
    k_idx =  logical(T.enabledAttribute(3,:));
    csc_names(~k_idx) = [];
    csc_names = csc_names+1; % offset o indexing.


    csc_list = dir([csc_dir filesep '*CH*.continuous']);
    for ii = length(csc_list):-1:1

        csc_name{ii} = csc_list(ii).name;
    end

    temp_list = [];

    for ii = 1:length(csc_idx)
        temp_idx = find(contains(csc_name, ['CH' num2str(csc_names(csc_idx(ii))) '.continuous']));
        if length(temp_idx) > 1
            error('Too many csc files containing this name')
        end
        if ~isempty(temp_idx)
            temp_list(ii).name= csc_name{temp_idx};
            temp_list(ii).folder= csc_list(1).folder;
        end
    end

    csc_list = temp_list;

    csc= []; labels = [];

    % get the first channel with the time vector.
    ii = 1;
    [data, tvec, info] = load_open_ephys_data([csc_list(ii).folder filesep csc_list(ii).name]);
    csc = tsd(tvec, data');
    labels{ii} = info.header.channel;
    csc.cfg.hdr{ii} = info.header;
    csc.cfg.hdr{ii}.SamplingFrequency = info.header.sampleRate;

    % hold the sampling frequency
    fs = csc.cfg.hdr{1}.SamplingFrequency;

    cfg_in.decimateFactor = 15;
    csc = decimate_tsd(cfg_in, csc);

    for ii = length(csc_list):-1:2
        [data, ~, info] = load_open_ephys_data([csc_list(ii).folder filesep csc_list(ii).name]);

        csc.data(ii,:) = decimate(data,cfg_in.decimateFactor)';

        labels{ii} = info.header.channel;
        csc.cfg.hdr{ii} = info.header;
        csc.cfg.hdr{ii}.SamplingFrequency = info.header.sampleRate ./ cfg_in.decimateFactor;
    end

    csc.label = labels;

    csc_prime = csc.tvec(1) + ts_prime;
    % csc.tvec = csc.tvec - csc.tvec(1);
    % csc.tvec = csc.tvec + ts_prime;
    csc.tvec = csc.tvec - csc.tvec(1) + (csc.tvec(1) - ts_prime(1)); % zero out the csc.
    % csc.tvec = csc.tvec + ts_prime(1); % zero out the csc.

    % cfg_in.decimateFactor = 15;
    % csc = decimate_tsd(cfg_in, csc);

    % load the OE version of the events.
    evts_list = dir([csc_dir filesep '*Data*.events']);

    OE_evts = OE_LoadEvents([evts_list.folder filesep evts_list.name], fs, ts_prime);


    % session info
    if ~contains(csc_dir, '/')
        parts = strsplit(csc_dir, '\');
    else
        parts = strsplit(csc_dir, '/');
    end
    s_idx = contains(parts, '2026'); % find the folder containing '2026' since all the sessions of interest are from this time.
    sess = parts{s_idx};
    if contains(sess, 'Pox')
        subject = lower(sess(strfind(sess, 'Pox'): strfind(sess, 'Pox')+6));
    elseif contains(lower(sess), 'som')
        subject = lower(sess(strfind(lower(sess), 'som'): strfind(lower(sess), 'som')+3));
    end
    if contains(sess, 'TFC')
        sess_id = strrep(sess(strfind(sess, 'TFC'):end), '_', '');
    else
        sess_id = sess(strfind(sess, 'LT'):end);
    end

    % make a quick test plot
    figure(10101)
    clf
    hold on
    if exist('S', 'var')
        plot(S)
    end
    plot(csc.tvec, csc.data(1,:)/100);

    if contains(OE_evts.label, '4')
        vline(OE_evts.t{ismember(OE_evts.label, '4')}, '--r')
    end
    if contains(OE_evts.label, '11')
        vline(OE_evts.t{ismember(OE_evts.label, '11')}, 'g')
    end
    if contains(OE_evts.label, '13')
        vline(OE_evts.t{ismember(OE_evts.label, '13')}, 'b')
    end
    if contains(OE_evts.label, '16')
        vline(OE_evts.t{ismember(OE_evts.label, '16')}, 'm')
    end
    vline(OE_evts.t{ismember(OE_evts.label, '8')}, 'c')

    xlim([csc.tvec(1)-20 csc.tvec(end)+20])



    %% make an average LFP channel if needed

    % csc.data(2,:) = mean(csc.data(2:end,:));
    % csc.label{2} = 'sub_mean';
    % csc.label(3:end) = [];
    % csc.data(3:end,:) = [];
    % csc.cfg.hdr(3:end) = [];
    %% csc check

    % restrict for speed

    csc_r = restrict(csc, csc.tvec(1), csc.tvec(1)+500);

    figure(100+iS)
    clf;
    hold on
    if exist('swr_idx', 'var')
        swr_ch = ismember(csc_r.label, swr_idx);
    else
        swr_ch = zeros(1,length(csc_r.label));
    end
    for ii = 1:size(csc_r.data,1)

        if swr_ch(ii)
            plot(csc_r.tvec, csc_r.data(ii,:)+ii*500,'k', 'LineWidth',1);
        else
            plot(csc_r.tvec, csc_r.data(ii,:)+ii*500);
        end

        lab{ii} = csc_r.label{ii};
        y_t(ii) = median(csc_r.data(ii,:)+ii*500);
    end
    set(gca, 'YTick', y_t, 'YTickLabel', lab)

        title(all_csc_dir{iS})

    % vline(SWR_evts.t{2}, 'r')

    % %% grab the SWR times
    % evts_list = dir([swr_dir filesep '*Data*.events']);
    % SWR_evts = OE_LoadEvents([evts_list.folder filesep evts_list.name], fs);
    %
    % iRi = diff(SWR_evts.t{2});
    % keep_idx = iRi <.05;
    %
    % SWR_evts.t{2}(keep_idx) = [];

    %% get the movement if present in the evts

    move_ts = ts({OE_evts.t{contains(OE_evts.label, '8')}}, {'move'});

    mov_rate  = MS_spike2rate(move_ts, csc.tvec);


    cfg_mov = [];
    cfg_mov.threshold = .001;
    cfg_mov.dcn = '<';
    cfg_mov.operation = '<';
    cfg_mov.minlen = .5;

    mov_iv = TSDtoIV(cfg_mov, mov_rate);

    % pad movement

    cfg_resize.d = [+.5 -.5];

    mov_iv = ResizeIV(cfg_resize, mov_iv);


    csc_r = restrict(csc, mov_iv);
    mov_r = restrict(mov_rate, mov_iv);

    figure(1000+iS)
    plot(csc_r.tvec, (csc_r.data(1,:)/100)-5)
    plot(mov_rate.tvec, ((mov_rate.data(1,:)/max(mov_rate.data(1,:)))*10)-10)
    ylim([-10 inf])
    title(all_csc_dir{iS})


    %% determine the probe type

    if sum(ismember(probe_type{2}, subject))>0
        probe{1} = '5x12'; probe{2} = 'Buz32'; 

        % reclassify spikes on shank 1 as being CA1 or DG
        
    elseif sum(ismember(probe_type{1}, subject))>0
        probe{1} = '4x16'; probe{2} = 'Buz32'; 
    end

    if strcmp(probe{1}, '5x12')
        dg_chan = (S.usr.ch > 8 & S.usr.ch < 17); 
        S_dg = SelectTS([], S, dg_chan);
        S_dg.loc = S.loc(dg_chan); 
        S_dg.c_ord = spring(sum(dg_chan));


        S = SelectTS([], S, ~dg_chan);
        S.loc = S.loc(~dg_chan); 
        S.c_ord(dg_chan,:) = [];
    end
    %% collect the data
    % load("all_TFC.mat")
    this_name = [subject '_' sess_id];

    all_TFC.(this_name).csc = csc;

    % all_TFC.(this_name).swrs_ca1= swrs_ca1;
    % all_TFC.(this_name).swrs_sub= swrs_sub;


    all_TFC.(this_name).mov_iv= mov_iv;
    all_TFC.(this_name).mov_vec= mov_rate;

    OE_evts.label = {'context' 'TFC_on' 'mov', 'tone1', 'tone2', 'lick', 'eye'};

    all_TFC.(this_name).evts=OE_evts;

    % all_TFC.(this_name).swr_idx = swr_idx;

    all_TFC.(this_name).phy_dir = phy_dir;
    all_TFC.(this_name).csc_dir = csc_dir;


    if exist('S', 'var')
        all_TFC.(this_name).S=S;
        all_TFC.(this_name).ts_prime = ts_prime;
        all_TFC.(this_name).csc_prime = csc_prime;
        all_TFC.(this_name).probe = probe;

        if exist('S_dg', 'var')
            all_TFC.(this_name).S_dg=S_dg;
        end
    end
    % save('all_TFC.mat', 'all_TFC')

    data = [];
    data.(this_name) = all_TFC.(this_name);

    save([this_name '.mat'], 'data', '-v7.3')

    disp('data saved')

    clearvars -except all_csc_dir iS data_dir inter_dir node probe_type csc_idx ca1_swr_idx to_cut
end






%% update with the SWRs



for  iS = 1
    %%   Detect CA1 SWRS

    swrs_ca1 = MS_SWR_detector(csc_r,swr_idx(1));

    % % manually select
    % swr_k = MS_manual_IV_selection(csc, swrs_ca1, .2);
    % swr_k_idx = ismember( 1:length(swrs_ca1.tend), swr_k)
    % swrs_ca1 = SelectIV([], swrs_ca1, logical(swr_k))
    %%   Detect Sub SWRS
    close all
    if length(swr_idx) > 2
        disp('Getting the mean lfp for Sub')

        this_idx = find(contains(csc_r.label, swr_idx));
        csc_r_m = csc_r;
        csc_r_m.data = mean(csc_r.data(this_idx,:));
        csc_r_m.label{1} = 'sub_mean';
        csc_r_m.label(2:end) = [];
        csc_r_m.cfg.hdr(2:end) = [];
        swrs_sub = MS_SWR_detector(csc_r_m,csc_r_m.label(1));
        %
    else

        swrs_sub = MS_SWR_detector(csc_r,swr_idx(2));
    end
    % swr_k = MS_manual_IV_selection(csc, swrs_sub, .2);
    % swr_k_idx = ismember( 1:length(swrs_sub.tend), swr_k)
    % swrs_sub = SelectIV([], swrs_sub, logical(swr_k_idx))

    %% get the deep vs superficial classifciation
    if contains(subject, 'SOM')
        probe{1} = '5x12';
        probe{2} = 'Buz32';
    elseif contains(subject, {'Pox2217','Pox3265', 'Pox3567', 'Pox3568'})
        probe{1} = 'A4x16';
        probe{2} = 'Buz32';
    end
    [S, rip_out] = HF_deep_super(S,csc, swrs_ca1, probe{1}, 1:64);
    S.usr.deep = double(S.usr.deep);
    S.usr.deep(~S.loc) = NaN;
    S.usr.rel_depth(~S.loc) = NaN;

end