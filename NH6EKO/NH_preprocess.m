function NH_preprocess(import_dir, export_dir)
%% NH_preprocess: load, append, and preprocess select channels for each subject and export to the inter_dir. Only needs to be run once.
%
%
%
%    Inputs:
%    - Requires the global PARAMS set up in Master_NH6EKO.m
%    - [optional] import_dir: [path]  path to raw data folder
%    - [optional] export_dir: [path]  path to export processed data
%
%
%    Outputs:
%    - none
%
%
%
%
% EC 2020-09-06   initial version
%
%
%
%% initialize

global PARAMS

if nargin == 0
    import_dir = PARAMS.raw_data_dir;
    export_dir = PARAMS.inter_dir;
end


%% loop through recording folders and load data for each subject.

% get the folders names in order
cd(import_dir)
this_dir = dir;
dirFlags = [this_dir.isdir];

this_dir = this_dir(dirFlags);
[~,idx] = sort([this_dir.datenum]);
this_dir = this_dir(idx);
keep_idx = zeros(1, length(this_dir));

for iD = 1:length(this_dir)
    if contains(this_dir(iD).name, 'NHE6KO')
        keep_idx(iD) = 1;
    end
end

this_dir = this_dir(logical(keep_idx));

this_dir.name

%% loop through for each subject

Subjects = fieldnames(PARAMS.Subjects);

for iSub = 1:length(Subjects)
    fprintf('<strong>%s</strong>: loading data for %s...\n', mfilename, Subjects{iSub});
    %     all_csc.tvec = [];
    %     all_csc.data = [];
    
    % only use the recordings blocks with the right dir_name for this Sub
    for iD = length(this_dir):-1:1
        if contains(this_dir(iD).name,PARAMS.Subjects.(Subjects{iSub}).dir_name)
            sub_keep_idx(iD) = 1;
        else
            sub_keep_idx(iD) = 0;
        end
    end
    this_sub_dir = this_dir(logical(sub_keep_idx));
    
    % run the recording blocks
    for iRec  = 1:length(this_sub_dir)
        cd(this_sub_dir(iRec).name) % move to the specific recording block dir.
        
        % load the csc for LFP and EMG channel decimate to 2k if needed.
        cfg_lfp.fc = {PARAMS.Subjects.(Subjects{iSub}).LFP_Chan}; % which channel to load
        cfg_lfp.desired_sampling_frequency = 2000; %desired sampling frequency.
        this_csc{iRec} = MS_LoadCSC(cfg_lfp);
        
        % load the emg
        if numel(PARAMS.Subjects.(Subjects{iSub}).EMG_Chan) ==1
            cfg_emg.fc = {PARAMS.Subjects.(Subjects{iSub}).EMG_Chan}; % which channel to load
        else
            cfg_emg.fc = PARAMS.Subjects.(Subjects{iSub}).EMG_Chan; % which channel to load
        end
        cfg_emg.desired_sampling_frequency = 2000; %desired sampling frequency.
        emg = MS_LoadCSC(cfg_emg);
        
        % determine how to process the emg.  
        
        % if infered through the xcorr of two LFP channels....
        if strcmp(PARAMS.Subjects.M11.emg_type, 'infered')
            emg_corr = MS_infer_EMG([], emg);
            
        else % or through the filtered raw channel; 
        % filter the EMG between 15 and 300Hz
        cfg_f_emg = [];
        cfg_f_emg.f = [15 300];
        cfg_f_emg.type = 'fdesign'; %the type of filter I want to use via filterlfp
        cfg_f_emg.order = 16; %type filter order
        emg_corr = FilterLFP(cfg_f_emg, emg_corr);
        
        end
        
        % put the emg back into the lfp struct as a data channel
        this_csc{iRec}.data(2,:) = emg.data;
        this_csc{iRec}.cfg.hdr{2} = emg.cfg.hdr{1};
        this_csc{iRec}.label = {'LFP', 'EMG'};
        
        clear emg
        % workaround for cheetah '5.6.3' which does not track the start time.  Had to get them from the log file.  
        if strcmp(Subjects{iSub}, 'M05') && iRec == 1
            this_csc{iRec}.cfg.hdr{1}.TimeCreated = '2020/09/04 144812';
        elseif strcmp(Subjects{iSub}, 'M05') && iRec == 2
            this_csc{iRec}.cfg.hdr{1}.TimeCreated = '2020/09/05 210032';
        end
        
        tstart = this_csc{iRec}.cfg.hdr{1}.TimeCreated;
        this_csc{iRec}.tstart = datetime(datestr([tstart(1:10) ' ' tstart(12:13) ':' tstart(14:15) ':' tstart(16:17)]));
        
        % add in events for recording start time vs acq start
        this_csc{iRec}.evt = LoadEvents([]);
        
        % append all the data
        %         all_csc.tvec = [all_csc.tvec; this_csc{iRec}.tvec];
        %         all_csc.data = [all_csc.data, this_csc{iRec}.data];
        cd(import_dir) % move back into the import dir.
    end % recording block
    
    
    % save each subject back as an intermediate
    fprintf('<strong>%s</strong>: saving data for %s...', mfilename, Subjects{iSub});
    save([PARAMS.inter_dir filesep Subjects{iSub} '_raw_data.mat'], 'this_csc', '-v7.3')
    %    save([PARAMS.inter_dir filesep Subjects{iSub} '_cat_data.mat'], 'all_csc', '-v7.3')
    
    clear all_csc this_csc
end% subject

