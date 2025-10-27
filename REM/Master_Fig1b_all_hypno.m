%% Master_REM_Sleep_stats

lfp_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Jisoo\Raw_LFP';
inter_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Hypno\LFP_data'; 
ca_data = [];

%% load the A_out to get the session information for all of the manuscript
% data

method = 'binary';

load(['C:\Users\ecarm\'  strrep('Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\inter\B_out_', '\', filesep) method '.mat'], 'B_out')

% rename due to progressive structure naming convention above.
A_out = B_out;
clear B_out

exclude_mouse = {'pv1254'};

for ii = length(A_out):-1:1

    if contains(A_out{ii}{1}.info.subject, exclude_mouse)
        fprintf('Removing sesson: <strong>%s</strong>\n', A_out{ii}{1}.info.subject);
        rm_idx(ii) = true;
    else
        rm_idx(ii) = false;
    end

end

A_out(rm_idx) = [];
A_list = [];
for iA = length(A_out):-1:1
    A_list{iA} = [upper(A_out{iA}{1}.info.subject) '_' A_out{iA}{1}.info.session];
end
%% collect the LFP data

L = dir(lfp_dir);
L_list = [];
for iL = length(L):-1:1
    if L(iL).isdir  && ~(strcmp(L(iL).name, '.') || strcmp(L(iL).name, '..'))
        L_list{iL} = L(iL).name;
    else
        L_list{iL} = [];
    end
end
L_list(cellfun(@isempty, L_list)) = []; % remove '.', '..' and non-dir;



for iA = 1:length(A_list)
    % ~contains(L_list, A_list{iA})
    idx = find(~cellfun(@isempty,strfind(L_list,A_list{iA}))); %
    if isempty(idx)
        disp(A_list{iA})
        continue
    else
        cd([lfp_dir filesep L_list{idx}]);

        % load the data and save an intermediate file
        cfg = [];
        cfg.desired_sampling_frequency = 2000;
        if contains(lower('PV1060'), lower(A_list{iA}(1:6)))
            cfg.fc = {'CSC1.ncs', 'CSC8.ncs'}; % CSC1 is the emg and the other is the best LFP channel
            emg_f = 0;
        elseif contains(lower('PV1252'), lower(A_list{iA}(1:6)))
            cfg.fc = {'CSC1.ncs', 'CSC6.ncs'}; % CSC1 is the emg and the other is the best LFP channel
            emg_f = 1;
        elseif contains(lower('PV1043'), lower(A_list{iA}(1:6)))
            cfg.fc = {'CSC1.ncs','CSC7.ncs'};
            emg_f = 0;
        else
            cfg.fc = {'CSC1.ncs','CSC7.ncs'};
            emg_f = 1;
        end

        out.evts = MS_LoadEvents();

        out.csc = MS_LoadCSC(cfg);

        % split the LFP and EMG
        lfp = out.csc;
        lfp.data = out.csc.data(2,:);
        lfp.label = out.csc.label{2};
        lfp.cfg.hdr= [];
        lfp.cfg.hdr{1} = out.csc.cfg.hdr{2};

        emg = out.csc;
        emg.data = out.csc.data(1,:);
        emg.label = out.csc.label{1};
        emg.cfg.hdr= [];
        emg.cfg.hdr{1} = out.csc.cfg.hdr{1};

        emg_thresh = 60;
        TD_ratio = 1.5;

        if emg_f
            cfg_emg = [];
            cfg_emg.threshold = 0;
            cfg_emg.f = [15 30];
            emg = FilterLFP(cfg_emg,emg);
        end
            out.hypno = MS_get_hypno(lfp, emg, [], emg_thresh, TD_ratio);

        out.hypno.emg_thresh = emg_thresh;
        out.hypno.td_r = TD_ratio;

        title(upper(A_list{iA}))
        save([inter_dir filesep 'Hypno_data_' upper(A_list{iA})],'out', '-mat')

        clearvars out csc emg lfp emg_data
    end

end

%% collect the data and make a summary plot

cd(inter_dir)

s_list = dir('*_data_*'); 

% loop over hypno sessions, restrict to the 2hours of sleep and save the
% hypnogram. 

for iS = length(s_list):-1:1

    load([inter_dir filesep s_list(iS).name])

    % clip data dim if not in TSD stardar nChan x nSample
    if size(out.hypno.data,1) > size(out.hypno.data,2)
        out.hypno.data = out.hypno.data'; 
    end

    % needed for tsd functions like 'restrict'
    out.hypno.type = 'tsd';
    out.hypno.cfg.history.mfun = []; 
    out.hypno.cfg.history.cfg =[]; 
    out.hypno.label = out.hypno.labels; % temporary needs to be fixed in main function. 

    out.hypno.data = downsample(out.hypno.data, 10); %decimate the data for speed and clarity from 2kHz to 200Hz; 
    out.hypno.tvec = downsample(out.hypno.tvec, 10); %decimate the data for speed and clarity from 2kHz to 200Hz; 

    % extract the hypnogram states and restrict
    start_rec_idx = contains(out.evts.label, 'Starting Recording'); 
    stop_rec_idx = contains(out.evts.label, 'Stopping Recording'); 

    pre_t = [out.evts.t{start_rec_idx}(1) out.evts.t{stop_rec_idx}(1)];
    post_t = [out.evts.t{start_rec_idx}(2) out.evts.t{stop_rec_idx}(2)];

    pre_hypno = restrict(out.hypno, pre_t(1), pre_t(2));
    post_hypno = restrict(out.hypno, post_t(1), post_t(2));

    all_pre.tvec{iS} = pre_hypno.tvec; 
    all_pre.data{iS} = pre_hypno.data; 
    all_post.tvec{iS} = post_hypno.tvec; 
    all_post.data{iS} = post_hypno.data; 
end

%% plot the data 

[max_pre, pre_idx]= max(cellfun(@length, all_pre.tvec)); 
[max_post,  post_idx]= max(cellfun(@length, all_pre.tvec)); 

pre_mat = zeros(length(all_pre.tvec), max_pre); 
post_mat = zeros(length(all_post.tvec), max_post); 

% populate a maxtrix with the hypnograms. 
for ii =1:length(all_pre.tvec)
    % session info
    p = strsplit(s_list(ii).name, '_');
    sess_info{ii} = [p{3}(3:end) '-' p{4}];

    l_pre = length(all_pre.tvec{ii}); 
    pre_mat(ii,end- l_pre+1:end) = all_pre.data{ii}; 

    l_post = length(all_post.tvec{ii}); 
    post_mat(ii,1:l_post) = all_post.data{ii}; 
end
pre_tvec = (all_pre.tvec{pre_idx} - all_pre.tvec{pre_idx}(end))/60 /60;
post_tvec = (all_post.tvec{post_idx} - all_post.tvec{post_idx}(1))/60 /60;

%
figure(1002)
clf
subplot(2,4,1)
imagesc(pre_tvec, 1:length(all_pre.tvec),pre_mat)
set(gca, 'YTick',1:length(all_pre.tvec), 'YTickLabel',sess_info)
colormap([1, 1, 1; viridis(3)]);  
cb = colorbar(gca, 'northoutside'); 
cb.Ticks = 0:3; cb.TickLabels = {'NA', 'wake', 'SWS', 'REM'};
ylabel('session')

subplot(2,4,2)
imagesc(post_tvec, 1:length(all_pre.tvec),post_mat)
set(gca, 'YTick',length(all_pre.tvec), 'YTickLabel',[])
cb = colorbar(gca, 'northoutside'); 
cb.Ticks = 0:3; cb.TickLabels = {'NA', 'wake', 'SWS', 'REM'}; 