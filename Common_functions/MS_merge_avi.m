function vid_out = MS_merge_avi(video_dir, save_dir, track_flag); 
%% MS_merge_avi:  loads the video files in a directory, merges them, applies tracking [optional] and saves an output. 
%
%
%
%    Inputs: 
%    -
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2024-10-25   initial version 
%
%
%
%% initialize

% get all the files in a directoy. 
cd(video_dir)
v_list = dir('*.avi');

% catch cases that do not start with a digit
matchC = regexp({v_list.name}, '^\d+\.avi', 'start');  % Or: '^[0-9]+\.txt'
match  = ~cellfun('isempty', matchC);
v_list   = v_list(match);

% sort the video list
v_num = []; 
for ii = 1:length(v_list)
v_num(ii) = str2double(v_list(ii).name(1:strfind(v_list(ii).name, '.avi')));

end
[~, s_idx] = sort(v_num); 

v_list = v_list(s_idx);




    j_files = dir('metaData.json');
    fid = fopen([j_files.name]);
    raw = fread(fid,inf);
    str = char(raw');
    meta_json = jsondecode(str);
    fclose(fid);
    
    meta_json.FS = str2double(meta_json.frameRate(1:2)); 
    
v = VideoWriter('MergedVideo.avi');% Create new video file
v.FrameRate =meta_json.FS*2; 
open(v)

for ii = 1:length(v_list)
    tic
    Vid1 = VideoReader(v_list(ii).name);
    
    while hasFrame(Vid1)
        frame = rgb2gray(readFrame(Vid1)); % read each frame
        writeVideo(v,frame) % write each frame
    end
        disp([v_list(ii).name])
    toc
end
close(v)