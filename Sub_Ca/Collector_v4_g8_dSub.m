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
    scope_id = 'My_v4_Miniscope';
    cam_id = 'My_WebCam';
elseif nargin <2
    inter_dir = [cd filesep 'Inter'];
    sub_id = 'dSub_g8';
    task_id = 'w_maze';
    scope_id = 'My_v4_Miniscope';
    cam_id = 'My_WebCam';
elseif nargin < 3
    sub_id = 'dSub_g8';
    task_id = 'w_maze';
    scope_id = 'My_v4_Miniscope';
    cam_id = 'My_WebCam';
elseif nargin <4
    task_id = 'w_maze';
    scope_id = 'My_v4_Miniscope';
    cam_id = 'My_WebCam';
elseif nargin <5
    scope_id = 'My_v4_Miniscope';
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
    fprintf('Subs found: %s\n', sub_list{ii})
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
        
        % check for processed data
        if ~exist([maze_dir filesep scope_id filesep 'minian_dataset.nc'], 'file') && ~exist([maze_dir filesep scope_id filesep 'minian_ms.mat'], 'file')
            fprintf('<strong> No Minian output found in</strong> %s', [maze_dir filesep scope_id])
            %            continue
        elseif exist([maze_dir filesep scope_id filesep 'minian_ms.mat'], 'file')
            fprintf('<strong> minian_ms found in</strong> %s', [maze_dir filesep scope_id])
            load([maze_dir filesep scope_id filesep 'minian_ms.mat'], 'ms');
            Maze.minianms = ms; clear ms; 
        elseif exist([maze_dir filesep scope_id filesep 'minian_dataset.nc'], 'file') && ~exist([maze_dir filesep scope_id filesep 'minian_ms.mat'], 'file')
            fprintf('<strong> minian_dataset.nc found in</strong> %s', [maze_dir filesep scope_id])
            
            Maze.minianms = Minian2MS([maze_dir filesep scope_id], 1, 1);
            saveas(gcf, [sub_id '_' sess_found{iS} '_Maze_ms.png'])
        end
        
        % check for processed HD
        if ~exist([maze_dir filesep scope_id filesep 'headOrientation.csv'], 'file')
            fprintf('<strong> No HD output found in</strong> %s', [maze_dir filesep scope_id])
            %            continue
        elseif exist([maze_dir filesep scope_id filesep 'headOrientation.csv'], 'file')
            fprintf('<strong> Processing HeadOrientation.csv found in</strong> %s', [maze_dir filesep scope_id])
            
            [Maze.HD, Maze.spd] = MS_Load_v4HD([maze_dir filesep scope_id], 0);
        end
        
         % check for processed tracking
        if ~exist([maze_dir filesep scope_id filesep '0_LocationOutput.csv'], 'file')
            fprintf('<strong> No tracking output found in</strong> %s', [maze_dir filesep scope_id])
            %            continue
        elseif exist([maze_dir filesep scope_id filesep '0_LocationOutput.csv'], 'file')
            fprintf('<strong> Processing 0_LocationOutput.csv found in</strong> %s', [maze_dir filesep scope_id])
            
            [Maze.pos, Maze.spd] = EZTrack2pos([maze_dir filesep scope_id], 1);
        end
        
        % interpolate to match data to ms. 
        if length(Maze.minianms.
        
        
    end
    
    
end
    
    
    
    
    
    
    
    
    
    
