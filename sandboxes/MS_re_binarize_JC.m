%% re-score cut data sandbox

% %% the raw traces are here for each recording block: 
% iRecord = 1; 
% ms_seg_resize.RawTraces{iRecord}

% Binarized was only save back into the individual segments.  but we can
% make new ones here. 

%Eric, before you make this script, i just wanto check, 
%so this ms_seg is the session that we already cut and it has raw traces
%that is cut so, i guess when we change the binarize, we can just run for
%each of this, right?
% that is correct.  It would be very simple to just have a loop that loads,
% reruns the binary with a new threshold and then saves it back.  The one
% thing that would need to happen is to redo all the 'all_binary_..."
% parts.  which would be easy to do as well. 
%Ok got it! i don't want to bother you for this, so i'll try to do it and
%upload my data! Or maybe guillaume might do it ? I just don't want to
%bother you haha

% writing the script would be very easy since it is just copying part of
% some other code and writing a quick loop.  One catch is that I can't do
% it tonight but I can have it done tomorrow.  it would not be any hassle. 
%Ok! if you have time then i would be appreciate it but i can also do it if
%you busy anyway. i'll just upload the data today ! and that's it!  Okay
%sounds good. 

%% get the file paths and loop

function ms_seg_resize = MS_re_binarize_JC(z_threshold, ms_load_dir, ms_save_dir, ms_fname_load, ms_fname_save)
%% MS_reclass_hynpo: simple script for reclassifying the hypno labels in segemented ms data.
%
%
%
%    Inputs:
%     - ms_load_dir [char]   path to ms_resize.mat file.
%
%     - ms_save_dir [char]   path to where you want the ms_resize.mat and the all_*.mat files to be updated and
%     saved.
%
%     - ms_fname_load [char] name of the .mat file to load.
%
%     - ms_fname_save [char] name of the .mat file to save
%

%    Outputs:
%     - ms_seg_resize: [struct] ms data that has been segmented.  contains
%     the updated ms_seg.resize.hypnolabel field.
%
%
%
%
% EC 2020-04-22   initial version
%
%
%% initialize

if nargin ==0
    z_threshold = 3; % default from ms_binarizeDetrend
    ms_load_dir = cd;
    ms_save_dir = cd;
    ms_fname_load = 'ms_resize'; % name of the file to load.
    ms_fname_save = 'ms_resize'; % name of the file to save.
    fprintf('<strong>%s</strong>: no dir inputs, using cd: %s \n', mfilename, cd);
elseif nargin == 1
    ms_load_dir = cd;
    ms_save_dir = ms_load_dir;
    ms_fname_load = 'ms_resize'; % name of the file to load.
    ms_fname_save = 'ms_resize'; % name of the file to save.
    fprintf('<strong>%s</strong>: only ms_load_dir specified, saving in same dir\n    ms_load_dir: %s\n    ms_save_dir: %s\n', mfilename, ms_load_dir, ms_load_dir);
elseif nargin == 2
    ms_save_dir = ms_load_dir;
    ms_fname_load = 'ms_resize'; % name of the file to load.
    ms_fname_save = 'ms_resize'; % name of the file to save.
    fprintf('<strong>%s</strong>: only ms_load_dir specified, saving in same dir\n    ms_load_dir: %s\n    ms_save_dir: %s\n', mfilename, ms_load_dir, ms_load_dir);
elseif nargin ==3 
    ms_fname_load = 'ms_resize'; % name of the file to load.
    ms_fname_save = 'ms_resize'; % name of the file to save.
    fprintf('<strong>%s</strong>: only ms_load_dir specified, saving in same dir\n    ms_load_dir: %s\n    ms_save_dir: %s\n', mfilename, ms_load_dir, ms_load_dir);
elseif nargin == 4
    ms_fname_save = 'ms_resize.mat'; % name of the file to save.
    fprintf('<strong>%s</strong>: only ms_load_dir specified, saving in same dir\n    ms_load_dir: %s\n    ms_save_dir: %s\n', mfilename, ms_load_dir, ms_save_dir);
end

%% load data
if exist([ms_fname_load '.mat'], 'file')
    load(ms_fname_load);
    
else
    error('no ms_resize.mat found')
end

%% load trk and update and extract time values. 
if exist(['ms_trk.mat'], 'file')
    load('ms_trk.mat');
    
else
    error('no ms_trk.mat found')
end

% get time
time_str = regexp(ms_trk.file_names,'\d*','Match');
ms_trk.time_labels = datestr(datestr([time_str{1} ':' time_str{2} ':' time_str{3}]), 'HH:MM:SS');
temp_time = datetime(ms_trk.time_labels); 
temp_time.Format = 'HH:mm:ss';
temp_time = temp_time+seconds((ms_trk.time(end)- ms_trk.time(1))/1000); 

trk_time = datevec(ms_trk.time_labels); 
trk_end_time = datevec(datestr(temp_time, 'HH:MM:SS')); 

% save the trk segment back with the updated time field. 
save('ms_trk.mat','ms_trk', '-v7.3')

% etime(trk_end_time, trk_time);
%% recompute the all pre/post SW vs REM blocks.

all_binary_pre = []; all_binary_post= [];
all_RawTraces_pre = []; all_RawTraces_post = [];
all_detrendRaw_pre = []; all_detrendRaw_post = [];
all_seg_idx = [];
%rem
all_binary_pre_REM = []; all_binary_post_REM= [];
all_RawTraces_pre_REM = []; all_RawTraces_post_REM = [];
all_detrendRaw_pre_REM = []; all_detrendRaw_post_REM = [];
pre_REM_idx = []; post_REM_idx = [];
% SW
all_binary_pre_SW = []; all_binary_post_SW= [];
all_RawTraces_pre_SW = []; all_RawTraces_post_SW = [];
all_detrendRaw_pre_SW = []; all_detrendRaw_post_SW = [];
pre_SW_idx = []; post_SW_idx = [];


for iSeg = 1:length(ms_seg_resize.RawTraces)
    ms_seg = []; % cleared so that we can use this var name for saving.
    
    keep_idx = 1:size(ms_seg_resize.RawTraces,1); % actually this is a remove index
    keep_idx =keep_idx(find((keep_idx ~= iSeg)));
    
    cfg_rem = [];
    ms_seg = MS_remove_data_sandbox(cfg_rem, ms_seg_resize, keep_idx);
    
    ms_seg = MS_de_cell(ms_seg);
    
    % binarize the trace
    
    ms_seg = MS_msExtractBinary_detrendTraces(ms_seg, z_threshold);
    
    % check for inactive cells and remove from ms.SFPs just using sum of
    % binary > 0;
    
    
    % add the spiking probability for different cell types.
    %     ms_seg.place_firing =
    
    % start of trk  =  trk_time; end of track = trk_end_time; 
    % get the time vector for current segment. 
    seg_time = datevec(ms_seg.time_labels); 
    % add time vector
    if strcmp(ms_seg.pre_post, 'pre')
         ms_seg.time2trk = etime(seg_time, trk_time)/60;
    elseif strcmp(ms_seg.pre_post, 'post')
         ms_seg.time2trk = etime(seg_time,trk_end_time)/60;
    end

    cfg_SFP = [];
    cfg_SFP.fnc = '==';
    cfg_SFP.remove_val = 0;
    ms_seg = MS_update_SFP(cfg_SFP, ms_seg);
    
    this_dir = [];
    this_dir = [ms_save_dir filesep ms_seg_resize.file_names{iSeg}];
    fprintf('<strong>%s</strong>: saving resized ms struct back to %s...\n', mfilename, this_dir)
    if exist(this_dir, 'dir') % if the dir already exists, delete it and all the subfiles/dir and make a new one.  Avoids problems with multiple ms_resize.mat files in one dir. 
        rmdir(this_dir, 's')
    end
    mkdir(this_dir)
    if  strcmp(ms_seg_resize.pre_post{iSeg}, 'post') && exist([this_dir filesep 'ms_seg_resize_' 'pre' '_' ms_seg_resize.hypno_label{iSeg} '.mat'])
        delete([this_dir filesep ms_fname_save '_' 'pre' '_' ms_seg_resize.hypno_label{iSeg} '.mat'])
    end
    % if this is a homecage do not use the 'pre' or post' lab.
    save([this_dir filesep ms_fname_save '_' ms_seg_resize.pre_post{iSeg} '_' ms_seg_resize.hypno_label{iSeg}],'ms_seg', '-v7.3');
    
    
    
    % keep the index for the segment.
    if isempty(all_seg_idx)
        all_seg_idx(iSeg) = length(ms_seg.RawTraces);
    else
        all_seg_idx(iSeg) = length(ms_seg.RawTraces) + all_seg_idx(iSeg -1);
    end
    % cat the binary traces for pre V post, and REM v SW
    if strcmp(ms_seg_resize.pre_post{iSeg}, 'pre')
        all_binary_pre = [all_binary_pre; ms_seg.Binary];
        all_RawTraces_pre = [all_RawTraces_pre; ms_seg.RawTraces];
        all_detrendRaw_pre = [all_detrendRaw_pre; ms_seg.detrendRaw];
        
        
        % break out REM and SW
        if strcmp(ms_seg_resize.hypno_label{iSeg}, 'REM')
            all_binary_pre_REM = [all_binary_pre_REM; ms_seg.Binary];
            all_RawTraces_pre_REM = [all_RawTraces_pre_REM; ms_seg.RawTraces];
            all_detrendRaw_pre_REM = [all_detrendRaw_pre_REM; ms_seg.detrendRaw];
            pre_REM_idx = [pre_REM_idx, iSeg];
        elseif strcmp(ms_seg_resize.hypno_label{iSeg}, 'SW')
            all_binary_pre_SW = [all_binary_pre_SW; ms_seg.Binary];
            all_RawTraces_pre_SW = [all_RawTraces_pre_SW; ms_seg.RawTraces];
            all_detrendRaw_pre_SW = [all_detrendRaw_pre_SW; ms_seg.detrendRaw];
            pre_SW_idx = [pre_SW_idx, iSeg];
        end
        
    elseif strcmp(ms_seg_resize.pre_post{iSeg}, 'post')
        all_binary_post = [all_binary_post; ms_seg.Binary];
        all_RawTraces_post = [all_RawTraces_post; ms_seg.RawTraces];
        all_detrendRaw_post = [all_detrendRaw_post; ms_seg.detrendRaw];
        
        % break out REM and SW
        if strcmp(ms_seg_resize.hypno_label{iSeg}, 'REM')
            all_binary_post_REM = [all_binary_post_REM; ms_seg.Binary];
            all_RawTraces_post_REM = [all_RawTraces_post_REM; ms_seg.RawTraces];
            all_detrendRaw_post_REM = [all_detrendRaw_post_REM; ms_seg.detrendRaw];
            post_REM_idx = [post_REM_idx, iSeg];
            
        elseif strcmp(ms_seg_resize.hypno_label{iSeg}, 'SW')
            all_binary_post_SW = [all_binary_post_SW; ms_seg.Binary];
            all_RawTraces_post_SW = [all_RawTraces_post_SW; ms_seg.RawTraces];
            all_detrendRaw_post_SW = [all_detrendRaw_post_SW; ms_seg.detrendRaw];
            post_SW_idx = [post_SW_idx, iSeg];
            
        end
    end
end
all_seg_idx = [0 all_seg_idx];


%% save the files. 
fprintf('<strong>%s</strong>: saving concatinating Binary, RawTraces, detrendRaw, and indicies\n', mfilename);

% save everything
save([ms_save_dir filesep 'all_seg_idx.mat'], 'all_seg_idx', '-v7.3');

% pre only
save([ms_save_dir filesep 'all_binary_pre.mat'], 'all_binary_pre', '-v7.3');
save([ms_save_dir filesep 'all_RawTraces_pre.mat'], 'all_RawTraces_pre', '-v7.3');
save([ms_save_dir filesep 'all_detrendRaw_pre.mat'], 'all_detrendRaw_pre', '-v7.3');

% post only
save([ms_save_dir filesep 'all_binary_post.mat' ], 'all_binary_post', '-v7.3');
save([ms_save_dir filesep 'all_RawTraces_post.mat'], 'all_RawTraces_post', '-v7.3');
save([ms_save_dir filesep 'all_detrendRaw_post.mat'], 'all_detrendRaw_post', '-v7.3');


% pre REM only
save([ms_save_dir filesep 'all_binary_pre_REM.mat'], 'all_binary_pre_REM', '-v7.3');
save([ms_save_dir filesep 'all_RawTraces_pre_REM.mat'], 'all_RawTraces_pre_REM', '-v7.3');
save([ms_save_dir filesep 'all_detrendRaw_pre_REM.mat'], 'all_detrendRaw_pre_REM', '-v7.3');

% post SW only
save([ms_save_dir filesep 'all_binary_post_REM.mat' ], 'all_binary_post_REM', '-v7.3');
save([ms_save_dir filesep 'all_RawTraces_post_REM.mat'], 'all_RawTraces_post_REM', '-v7.3');
save([ms_save_dir filesep 'all_detrendRaw_post_REM.mat'], 'all_detrendRaw_post_REM', '-v7.3');

% pre REM only
save([ms_save_dir filesep 'all_binary_pre_SW.mat'], 'all_binary_pre_SW', '-v7.3');
save([ms_save_dir filesep 'all_RawTraces_pre_SW.mat'], 'all_RawTraces_pre_SW', '-v7.3');
save([ms_save_dir filesep 'all_detrendRaw_pre_SW.mat'], 'all_detrendRaw_pre_SW', '-v7.3');

% post SW only
save([ms_save_dir filesep 'all_binary_post_SW.mat' ], 'all_binary_post_SW', '-v7.3');
save([ms_save_dir filesep 'all_RawTraces_post_SW.mat'], 'all_RawTraces_post_SW', '-v7.3');
save([ms_save_dir filesep 'all_detrendRaw_post_SW.mat'], 'all_detrendRaw_post_SW', '-v7.3');

%% save a file with the zscore config
% save([ms_save_dir filesep ms_fname_save '.mat'], 'ms_seg_resize', '-v7.3')
cfgs = [];
cfgs.z_theshold = z_threshold; 
cfgs.date = datestr(date, 'yyyy_mm_dd'); 
save([ms_save_dir filesep 'cfgs_z_' strrep(num2str(z_threshold), '.','p') '_on_' cfgs.date   '.mat'], 'cfgs', '-v7.3')


end % end function