function ms_seg_resize = MS_homecage_class(first_post_fname, ms_seg_resize)
%% MS_homecage_class: Reclassifies the output from MS_Segment_raw (_EV)
%  into pre and post blocks based on user input.  To run move to the
%  resized ms directory and run this script either with the block that
%  begins the 'post' period as an input or leave it empty to have the
%  script list the blocks and prompt the user for input.
%
%
%
%    Inputs:
%     - first_post_fname: [string] the name of the first block of the
%     'post' period. Should match the original ms data block folder name
%     (ie: 'H10_M20_S42_SWS60s')
%
%   OPTIONAL
%    - ms_seg_resize: [struct] the segmented ms struct from MS_Segment_raw.
%     Uses a specific ms struct rather than just loading the one in the
%     current dir.
%
%
%    Outputs:
%     - ms_seg_resize_out: updated ms struct with the correct pre and post
%     labels.
%
%
%
%
% EC 2020-04-13   initial version
%
%
%%  initialize
% user can specify the first_post_name and the ms struct if they wish.
% Otherwise it will load the local ms_resize and prompt the user for the
% file name to use when delineating.


if nargin == 0
    user_prompt = 1;
    ms_seg_resize = [];
    
    fprintf('<strong>%s</strong>: Using user prompt to delineate ''pre'' and ''post'' in local ms_resize.mat\n', mfilename)
    if ~exist('ms_resize.mat')
        error('No ms_resize.mat found in this dir')
    end
elseif nargin == 1
    user_prompt = 0;
    ms_seg_resize = [];
    
    fprintf('<strong>%s</strong>: Using <strong>%s</strong> to delineate ''pre'' and ''post'' in local ms_resize.mat\n', mfilename, first_post_fname)
    if ~exist('ms_resize.mat')
        error('No ms_resize.mat found in this dir')
    end
else
    user_prompt = 0;
    fprintf('<strong>%s</strong>: Using <strong>%s</strong> to delineate ''pre'' and ''post'' in user input ms_seg_resize\n', mfilename, first_post_fname)
end


%% Get the ms_resize.mat

if isempty(ms_seg_resize)
    load('ms_resize.mat');
end

%% get the first_post_fname: this is the name of the file
if user_prompt ==1
    fprintf('<strong>%s</strong> file names: \n', 'ms_seg_resize')
    for iB = 1:length(ms_seg_resize.file_names)
        fprintf(' %d)     <strong>%s</strong>\n', iB, ms_seg_resize.file_names{iB})
    end
    
    reply = input('\n\n<strong>Which file name is the first of the ''post'' period? use either full name or the index is the above list</strong>\n\n');
    
    if ischar(reply)
        first_post_fname = reply;
        fprintf('\n\nUsing <strong>%s</strong> as the first ''post'' of the post period\n', mfilename, first_post_fname)
    elseif isnumeric(reply)
        first_post_fname = ms_seg_resize.file_names{reply};
        fprintf('\n\n[numeric input] Using <strong>%s</strong> as the first ''post'' of the post period\n', first_post_fname)
    end
end


%% classify
% which idx is the start of the ''post' period
post_start_idx = find(strcmp(first_post_fname, ms_seg_resize.file_names));

for iB = 1:length(ms_seg_resize.file_names)
    
    if iB < post_start_idx
        fprintf(' %d)  <strong>%s</strong>  %s\n', iB, 'pre', ms_seg_resize.file_names{iB})
        ms_seg_resize.pre_post{iB} = 'pre';
    elseif iB == post_start_idx
        fprintf('____________________________________\n')
        fprintf(' %d)        <strong>%s  %s</strong>\n', iB, 'post', ms_seg_resize.file_names{iB})
        ms_seg_resize.pre_post{iB} = 'post';
    else
        fprintf(' %d)       <strong>%s</strong>  %s\n', iB, 'post', ms_seg_resize.file_names{iB})
        ms_seg_resize.pre_post{iB} = 'post';
    end
end



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
    
    ms_seg = msExtractBinary_detrendTraces(ms_seg);
    
    % check for inactive cells and remove from ms.SFPs just using sum of
    % binary > 0;
    
    cfg_SFP = [];
    cfg_SFP.fnc = '==';
    cfg_SFP.remove_val = 0;
    ms_seg = MS_update_SFP(cfg_SFP, ms_seg);
    
    this_dir = [];
    this_dir = [cd filesep ms_seg_resize.file_names{iSeg}];
    fprintf('<strong>%s</strong>: saving resized ms struct back to %s...\n', mfilename, this_dir)
    mkdir(this_dir)
    if  strcmp(ms_seg_resize.pre_post{iSeg}, 'post') && exist([this_dir filesep 'ms_seg_resize_' 'pre' '_' ms_seg_resize.hypno_label{iSeg} '.mat'])
        delete([this_dir filesep 'ms_seg_resize_' 'pre' '_' ms_seg_resize.hypno_label{iSeg} '.mat'])
    end
    % if this is a homecage do not use the 'pre' or post' lab.
    save([this_dir filesep 'ms_seg_resize_' ms_seg_resize.pre_post{iSeg} '_' ms_seg_resize.hypno_label{iSeg}],'ms_seg', '-v7.3');
    
    
    
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

fprintf('<strong>%s</strong>: saving concatinating Binary, RawTraces, detrendRaw, and indicies\n', mfilename);

% save everything
save([cd filesep 'all_seg_idx.mat'], 'all_seg_idx', '-v7.3');

% pre only
save([cd filesep 'all_binary_pre.mat'], 'all_binary_pre', '-v7.3');
save([cd filesep 'all_RawTraces_pre.mat'], 'all_RawTraces_pre', '-v7.3');
save([cd filesep 'all_detrendRaw_pre.mat'], 'all_detrendRaw_pre', '-v7.3');

% post only
save([cd filesep 'all_binary_post.mat' ], 'all_binary_post', '-v7.3');
save([cd filesep 'all_RawTraces_post.mat'], 'all_RawTraces_post', '-v7.3');
save([cd filesep 'all_detrendRaw_post.mat'], 'all_detrendRaw_post', '-v7.3');


% pre REM only
save([cd filesep 'all_binary_pre_REM.mat'], 'all_binary_pre_REM', '-v7.3');
save([cd filesep 'all_RawTraces_pre_REM.mat'], 'all_RawTraces_pre_REM', '-v7.3');
save([cd filesep 'all_detrendRaw_pre_REM.mat'], 'all_detrendRaw_pre_REM', '-v7.3');

% post SW only
save([cd filesep 'all_binary_post_REM.mat' ], 'all_binary_post_REM', '-v7.3');
save([cd filesep 'all_RawTraces_post_REM.mat'], 'all_RawTraces_post_REM', '-v7.3');
save([cd filesep 'all_detrendRaw_post_REM.mat'], 'all_detrendRaw_post_REM', '-v7.3');

% pre REM only
save([cd filesep 'all_binary_pre_SW.mat'], 'all_binary_pre_SW', '-v7.3');
save([cd filesep 'all_RawTraces_pre_SW.mat'], 'all_RawTraces_pre_SW', '-v7.3');
save([cd filesep 'all_detrendRaw_pre_SW.mat'], 'all_detrendRaw_pre_SW', '-v7.3');

% post SW only
save([cd filesep 'all_binary_post_SW.mat' ], 'all_binary_post_SW', '-v7.3');
save([cd filesep 'all_RawTraces_post_SW.mat'], 'all_RawTraces_post_SW', '-v7.3');
save([cd filesep 'all_detrendRaw_post_SW.mat'], 'all_detrendRaw_post_SW', '-v7.3');

%% save the resized ms back in with the updated pre/post field. 
save([cd filesep 'ms_resize.mat'], 'ms_seg_resize', '-v7.3')


%% visualize
figure(1010)
subplot(3,1,1)
title('Binary Pre cat: SW: green REM: red')
hold on
for ii = 1:10
    plot(all_binary_pre(:,ii)+ii);
end
xlim([0 length(all_binary_pre)])
% put in vertical lines for the start of REm or SW blocks, if they are
% present. 
% if ~isempty(pre_REM_idx)
%     vline(all_seg_idx(pre_REM_idx), {'r'});
% end
% if ~isempty(pre_SW_idx)
%     vline(all_seg_idx(pre_SW_idx),{'g'});
% end
subplot(3,1,2)
title('RawTraces Pre cat: SW: green REM: red')
hold on
for ii = 1:10
    plot(all_RawTraces_pre(:,ii)+ii);
end
% vline(all_seg_idx(pre_REM_idx), {'r'});
% vline(all_seg_idx(pre_SW_idx),{'g'})
xlim([0 length(all_RawTraces_pre)])

subplot(3,1,3)
title('detrendRaw Pre cat: SW: green REM: red')

hold on
for ii = 1:10
    plot(all_detrendRaw_pre(:,ii)+ii);
end
% vline(all_seg_idx(pre_REM_idx), {'r'});    
% vline(all_seg_idx(pre_SW_idx),{'g'})
xlim([0 length(all_detrendRaw_pre)])

pause(5)

saveas(gcf, [cd filesep 'cat_check'],'fig')
saveas(gcf, [cd filesep 'cat_check'],'png')

close; 
fprintf('<strong>%s</strong>: complete and saved\n', mfilename)
