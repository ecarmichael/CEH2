function data = sandbox_TFC_CA_collect(sess_dir)


cd(sess_dir)

%% loop over recording blocks

mlist = dir(fullfile(cd, ['**' filesep 'SFP_k.mat']));  %get list of files and folders in any subfolder
mlist = mlist(~[mlist.isdir]);


% loop and find the order of the recordings;
r_t = [];
for ii = 1:length(mlist)

    this_json = MS_load_json([MS_parent_dir(mlist(ii).folder) filesep 'metaData.json']);

    % disp([mlist(ii).folder ' - ' this_json.recordingStartTime   ])
    r_t(ii) = this_json.recordingStartTime.msec;
end

% sort based on the internal clock at recording start.
[~, f_s] = sort(r_t, "ascend");

mlist = mlist(f_s);

%%
for ii = 1:length(mlist)

    disp([mlist(ii).folder filesep mlist(ii).name])

    % load the ms (should have keep_idx)
    warning off
    load([mlist(ii).folder filesep 'ms.mat'])
    warning on

    % if it has been currated then remove the bad cell.
    if ~isfield(ms, 'keep_idx')
        fprintf('%s does not have a keep_idx', [mlist(ii).folder filesep 'ms'])
    else
        ms = MS_Ca_good_cells(ms);
    end

    % reduce size of ms
    ms_r = ms;
    ms_r = rmfield(ms, {'numFiles', 'numFrames', 'vidNum', 'frameNum', 'maxFramesPerFile', 'vidObj', 'dateNum', 'PeakToNoiseProj','analysis_time'});

    ms = ms_r; clear ms_r;

    % append the json
    sess_json = MS_load_json([MS_parent_dir(mlist(ii).folder) filesep 'metaData.json']);
    scope_json = MS_load_json([mlist(ii).folder filesep 'metaData.json']);

    % get the hd
    hd = MS_Load_v4HD(mlist(ii).folder);

    if ii == 1
        data.pre.ms = ms;
        data.pre.hd = hd;
        args=[namedargs2cell(sess_json), namedargs2cell(scope_json)];
        data.pre.json = struct(args{:});
    elseif ii == 2
        data.task.ms = ms;
        data.task.hd = hd;
        args=[namedargs2cell(sess_json), namedargs2cell(scope_json)];
        data.task.json = struct(args{:});
    elseif ii == 3
        data.post.ms = ms;
        data.post.hd = hd;
        args=[namedargs2cell(sess_json), namedargs2cell(scope_json)];
        data.post.json = struct(args{:});
    end
    clear ms hd sess_json scope_json

end

