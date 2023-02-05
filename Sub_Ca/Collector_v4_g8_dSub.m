function Collector_v4_g8_dSub(data_dir, inter_dir,sub_id, task_id, scope_id, cam_id)
%% Collector_v4_g8_dSub: Preprocess and collect v4 data for the dSub g8 W_maze and OF data.
%  - data structs:
%       - HD from 'head
%
%
%    Inputs:
%    - data_dir: [path]  path to data folders
%
%    - inter_dir: [path] where you want the intermediate files to be saved.
%
%    - sub_id: [string]  common string to all subjects to include. Ex:
%    'dSub_g8'
%    - task_id: [string]  name of the task or experiment. default is
%    'w_maze'.
%
%    - scope_id: [string] name of the miniscope. default is
%    'My_v4_Miniscope'
%
%    - cam_id: [string] name of the behav cam. default is
%    'My_webCam'
%
%    Outputs:
%    -
%
%
%
%
% EC 2023-01-25   initial version
%
%
%
%% initialize

if nargin <1
    data_dir = cd;
    inter_dir = [cd filesep 'Inter'];
    sub_id = 'dSub_g8';
    task_id = 'w_maze';
    scope_id = 'My_V4_Miniscope';
    cam_id = 'My_WebCam';
elseif nargin <2
    inter_dir = [cd filesep 'Inter'];
    sub_id = 'dSub_g8';
    task_id = 'w_maze';
    scope_id = 'My_V4_Miniscope';
    cam_id = 'My_WebCam';
elseif nargin < 3
    sub_id = 'dSub_g8';
    task_id = 'w_maze';
    scope_id = 'My_V4_Miniscope';
    cam_id = 'My_WebCam';
elseif nargin <4
    task_id = 'w_maze';
    scope_id = 'My_V4_Miniscope';
    cam_id = 'My_WebCam';
elseif nargin <5
    scope_id = 'My_V4_Miniscope';
    cam_id = 'My_WebCam';
elseif nargin <6
    cam_id = 'My_WebCam';
    
    
end

cd(data_dir);
fprintf('Looking %s %s %s in: %s ...\n',sub_id, task_id, scope_id, data_dir)
fprintf('Saving  in: %s ...\n', inter_dir)



%% List subjects matching the criteria


subs_found = dir([sub_id '*']);

sub_list = [];
for ii = length(subs_found):-1:1
    sub_list{ii} = subs_found(ii).name;
    fprintf('Subjects found: <strong>%s</strong>\n', sub_list{ii})
end


% loop over subjects and find all sessions.

for iSub = 1:length(sub_list)
    this_sub_dir = [data_dir filesep sub_list{iSub} filesep task_id];
    cd(this_sub_dir)
    sess_found = dir('2*');
    
    for iS = 1:length(sess_found)
        cd([this_sub_dir filesep sess_found(iS).name])
        recs_found = dir;
        % remove the first two which are '.' and '..'.
        recs_list = {recs_found([recs_found.isdir]).name};
        recs_list = recs_list(~ismember(recs_list ,{'.','..'}))';
        
        % open the timestamps from each recording block and sort them by
        % duration.
        rec_dir = []; rec_len = [];
        for iR = length(recs_list):-1:1
            rec_dir{iR} = [this_sub_dir filesep sess_found(iS).name filesep recs_list{iR}];
            this_TS = readtable([rec_dir{iR} filesep scope_id filesep 'timeStamps.csv']);
            
            this_TS = table2array(this_TS(:,2));
            
            rec_len(iR) = ((this_TS(end) - this_TS(1))/1000)/60;
        end
        
        [rec_len, sort_idx] = sort(rec_len, 'descend');
        
        rec_dir = rec_dir(sort_idx);
        
        maze_dir = rec_dir{1};
        OF_dir = rec_dir{2};
        
        fprintf('\n<strong>W maze %.2fmin dir:</strong> %s\n',rec_len(1), maze_dir)
        fprintf('<strong>OF     %.2fmin dir:</strong> %s\n',rec_len(2), OF_dir)
        
        %% Maze data
        % check for processed data
        if ~exist([maze_dir filesep scope_id filesep 'minian_dataset.nc'], 'file') && ~exist([maze_dir filesep scope_id filesep 'minian_ms.mat'], 'file')
            fprintf('<strong> No Minian output found in</strong> %s', [maze_dir filesep scope_id])
%             continue
        elseif exist([maze_dir filesep scope_id filesep 'minian_ms.mat'], 'file')
            fprintf('<strong> minian_ms found in</strong> %s', [maze_dir filesep scope_id])
            load([maze_dir filesep scope_id filesep 'minian_ms.mat'], 'ms');
            Maze.minianms = ms; clear ms;
        elseif exist([maze_dir filesep scope_id filesep 'minian_dataset.nc'], 'file') && ~exist([maze_dir filesep scope_id filesep 'minian_ms.mat'], 'file')
            fprintf('<strong> minian_dataset.nc found in</strong> %s', [maze_dir filesep scope_id])
            
            Maze.minianms = Minian2MS([maze_dir filesep scope_id], 1, 1);
            %             saveas(gcf, [sub_id '_' sess_found(iS).name '_Maze_ms.png'])
        end
        
        % check for processed HD
        if ~exist([maze_dir filesep scope_id filesep 'headOrientation.csv'], 'file')
            fprintf('<strong> No HD output found in</strong> %s', [maze_dir filesep scope_id])
            %            continue
        elseif exist([maze_dir filesep scope_id filesep 'headOrientation.csv'], 'file')
            fprintf('<strong> Processing HeadOrientation.csv found in</strong> %s', [maze_dir filesep scope_id])
            
            [Maze.HD] = MS_Load_v4HD([maze_dir filesep scope_id], 0);
        end
        
        % check for processed tracking
        if isempty(dir([maze_dir filesep cam_id filesep '*DLC*.csv']))
            fprintf('<strong> No DLC tracking output found in</strong> %s', [maze_dir filesep cam_id])
            %            continue
        elseif ~isempty(dir([maze_dir filesep cam_id filesep '*DLC*.csv']))
            DLC_fname = dir([maze_dir filesep cam_id filesep '*DLC*.csv']);
            fprintf('<strong> Processing %s found in</strong> %s',DLC_fname(1).name, [maze_dir filesep cam_id])
            
            Maze.pos = MS_DLC2TSD([maze_dir filesep cam_id], [], 0);
            %             [Maze.pos, Maze.spd] = EZTrack2pos([maze_dir filesep scope_id], 1);
            
            % fill in points outside of ROI
            figure(1010)
            clf
            maximize
            plot(Maze.pos.data(1,:), Maze.pos.data(2,:), '.');
            maze_out_roi = drawpolygon(gca, 'color', 'blue');
            
            maze_in_roi = drawpolygon(gca, 'color', 'red');
            maze_in2_roi = drawpolygon(gca, 'color', 'red');
            
            in_idx = inpolygon(Maze.pos.data(1,:), Maze.pos.data(2,:), maze_out_roi.Position(:,1), maze_out_roi.Position(:,2));
            out_idx = inpolygon(Maze.pos.data(1,:), Maze.pos.data(2,:), maze_in_roi.Position(:,1), maze_in_roi.Position(:,2));
            if exist('maze_in2_roi')
                out2_idx = inpolygon(Maze.pos.data(1,:), Maze.pos.data(2,:), maze_in2_roi.Position(:,1), maze_in2_roi.Position(:,2));
                out_idx = out_idx | out2_idx;
            end
            
            % check
            figure(1011)
            clf
            hold on
            %             plot(Maze.pos.data(1,:), Maze.pos.data(2,:), 'b.');
            plot(Maze.pos.data(1,in_idx & ~out_idx),Maze.pos.data(2,in_idx & ~out_idx),'b.') % points inside
            plot(Maze.pos.data(1,~in_idx | out_idx),Maze.pos.data(2,~in_idx | out_idx),'rx') % points outside
            
            % remove the points outside the boundries and replace with
            % nearest neighbour.
            Maze.pos.data(1:2,~in_idx | out_idx) = NaN;
            Maze.pos.data(1,:) = fillmissing(Maze.pos.data(1,:), 'nearest');
            Maze.pos.data(2,:) = fillmissing(Maze.pos.data(2,:), 'nearest');
            plot(Maze.pos.data(1,:), Maze.pos.data(2,:), 'bo');
            legend({'Inbounds', 'Outbounds', 'Missing Filled'});
            pause(1)
            close(1011)
            close(1010)
        end
        
        % interpolate to match data to ms.
        
        ms_length = length(Maze.minianms.RawTraces);
        ms_times = [Maze.minianms.tvec(1) Maze.minianms.tvec(end)];
        
        hd_length = length(Maze.HD.data);
        hd_times = [Maze.HD.tvec(1) Maze.HD.tvec(end)];
        
        if hd_length ~= ms_length
            fprintf('Length of data differs between <strong>%s (%.3f - %.3fs)</strong> and <strong>%s (%.3f - %.3fs)</strong>\n', 'Minian_ms', ms_times(1), ms_times(2), 'HD', hd_times(1), hd_times(2))
            int_data = [];
            for ii = size(Maze.HD.data,1):-1:1
                int_data(ii,:) = interp1(Maze.HD.tvec,Maze.HD.data(ii,:),  Maze.minianms.tvec);
            end
            Maze.HD.data = int_data;
        end
        Maze.HD.tvec = Maze.minianms.tvec;
        
        % position data.
        pos_length = length(Maze.pos.data);
        pos_times = [Maze.pos.tvec(1) Maze.pos.tvec(end)];
        if pos_length ~= ms_length
            fprintf('Length of data differs between <strong>%s (%.3f - %.3fs)</strong> and <strong>%s (%.3f - %.3fs)</strong>\n', 'Minian_ms', ms_times(1), ms_times(2), 'Pos', pos_times(1), pos_times(2))
            int_data = [];
            for ii = size(Maze.pos.data,1):-1:1
                int_data(ii,:) = interp1( Maze.pos.tvec, Maze.pos.data(ii,:), Maze.minianms.tvec);
            end
            Maze.pos.data = int_data;
        end
        Maze.pos.tvec = Maze.minianms.tvec;
        
        
        %% Open Field
        
        % check for processed data
        if ~exist([OF_dir filesep scope_id filesep 'minian_dataset.nc'], 'file') && ~exist([OF_dir filesep scope_id filesep 'minian_ms.mat'], 'file')
            fprintf('<strong> No Minian output found in</strong> %s', [OF_dir filesep scope_id])
            %            continue
        elseif exist([OF_dir filesep scope_id filesep 'minian_ms.mat'], 'file')
            fprintf('<strong> minian_ms found in</strong> %s', [OF_dir filesep scope_id])
            load([OF_dir filesep scope_id filesep 'minian_ms.mat'], 'ms');
            OF.minianms = ms; clear ms;
        elseif exist([OF_dir filesep scope_id filesep 'minian_dataset.nc'], 'file') && ~exist([OF_dir filesep scope_id filesep 'minian_ms.mat'], 'file')
            fprintf('<strong> minian_dataset.nc found in</strong> %s', [OF_dir filesep scope_id])
            
            OF.minianms = Minian2MS([OF_dir filesep scope_id], 1, 1);
            %             saveas(gcf, [sub_id '_' sess_found(iS).name '_OF_ms.png'])
        end
        
        % check for processed HD
        if ~exist([OF_dir filesep scope_id filesep 'headOrientation.csv'], 'file')
            fprintf('<strong> No HD output found in</strong> %s', [OF_dir filesep scope_id])
            %            continue
        elseif exist([OF_dir filesep scope_id filesep 'headOrientation.csv'], 'file')
            fprintf('<strong> Processing HeadOrientation.csv found in</strong> %s', [OF_dir filesep scope_id])
            
            [OF.HD] = MS_Load_v4HD([OF_dir filesep scope_id], 0);
        end
        
        % check for processed tracking
        if isempty(dir([OF_dir filesep cam_id filesep '*DLC*.csv']))
            fprintf('<strong> No DLC tracking output found in</strong> %s', [OF_dir filesep cam_id])
            %            continue
        elseif ~isempty(dir([OF_dir filesep cam_id filesep '*DLC*.csv']))
            DLC_fname = dir([OF_dir filesep cam_id filesep '*DLC*.csv']);
            fprintf('<strong> Processing %s found in</strong> %s',DLC_fname(1).name, [OF_dir filesep cam_id])
            
            OF.pos = MS_DLC2TSD([OF_dir filesep cam_id], [], 0);
            %             [OF.pos, OF.spd] = EZTrack2pos([OF_dir filesep scope_id], 1);
            
            % fill in points outside of ROI
            figure(1010)
            clf
            maximize
            plot(OF.pos.data(1,:), OF.pos.data(2,:), '.');
            OF_out_roi = drawpolygon(gca, 'color', 'blue');
            
            
            in_idx = inpolygon(OF.pos.data(1,:), OF.pos.data(2,:), OF_out_roi.Position(:,1), OF_out_roi.Position(:,2));
            
            
            % check
            figure(1011)
            clf
            hold on
            %             plot(OF.pos.data(1,:), OF.pos.data(2,:), 'b.');
            plot(OF.pos.data(1,in_idx ),OF.pos.data(2,in_idx),'b.') % points inside
            plot(OF.pos.data(1,~in_idx),OF.pos.data(2,~in_idx ),'rx') % points outside
            
            % remove the points outside the boundries and replace with
            % nearest neighbour.
            OF.pos.data(1:2,~in_idx ) = NaN;
            OF.pos.data(1,:) = fillmissing(OF.pos.data(1,:), 'nearest');
            OF.pos.data(2,:) = fillmissing(OF.pos.data(2,:), 'nearest');
            plot(OF.pos.data(1,:), OF.pos.data(2,:), 'bo');
            legend({'Inbounds', 'Outbounds', 'Missing Filled'});
                 pause(1)
            close(1011)
            close(1010)
        end
        
        if exist('OF', 'var')
            if sum(contains(fieldnames(OF), {'minianms','HD', 'pos'})) ~=3
                fprintf(['<strong>' sub_list{iSub} '_' sess_found(iS).name ' is incomplete</strong>\n'])
                continue
            end
        end
        % interpolate to match data to ms.
        
        ms_length = length(OF.minianms.RawTraces);
        ms_times = [OF.minianms.tvec(1) OF.minianms.tvec(end)];
        
        hd_length = length(OF.HD.data);
        hd_times = [OF.HD.tvec(1) OF.HD.tvec(end)];
        
        if hd_length ~= ms_length
            fprintf('Length of data differs between <strong>%s (%.3f - %.3fs)</strong> and <strong>%s (%.3f - %.3fs)</strong>\n', 'Minian_ms', ms_times(1), ms_times(2), 'HD', hd_times(1), hd_times(2))
            int_data = [];
            for ii = size(OF.HD.data,1):-1:1
                int_data(ii,:) = interp1(OF.HD.tvec,OF.HD.data(ii,:),  OF.minianms.tvec);
            end
            OF.HD.data = int_data;
        end
        OF.HD.tvec = OF.minianms.tvec;
        
        % position data.
        pos_length = length(OF.pos.data);
        pos_times = [OF.pos.tvec(1) OF.pos.tvec(end)];
        if pos_length ~= ms_length
            fprintf('Length of data differs between <strong>%s (%.3f - %.3fs)</strong> and <strong>%s (%.3f - %.3fs)</strong>\n', 'Minian_ms', ms_times(1), ms_times(2), 'Pos', pos_times(1), pos_times(2))
            int_data = [];
            for ii = size(OF.pos.data,1):-1:1
                int_data(ii,:) = interp1( OF.pos.tvec, OF.pos.data(ii,:), OF.minianms.tvec);
            end
            OF.pos.data = int_data;
        end
        OF.pos.tvec = OF.minianms.tvec;
        
        %% save the data back to the intermediate dir
        if ~exist(inter_dir, 'dir')
            mkdir(inter_dir)
        end
        save([inter_dir filesep sub_list{iSub} '_' sess_found(iS).name '_maze.mat'],'Maze');
        save([inter_dir filesep sub_list{iSub} '_' sess_found(iS).name '_OF.mat'],'OF');
        
        clear Maze OF
    end
    
    
end










