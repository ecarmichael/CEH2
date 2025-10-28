%% Master_REM_Sleep_stats

lfp_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Jisoo\Raw_LFP';

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
            emg_f = 1;
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
            cfg_emg.f = [10 20];
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