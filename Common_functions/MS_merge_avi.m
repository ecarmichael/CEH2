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


    j_files = dir('metaData.json');
    fid = fopen([j_files.name]);
    raw = fread(fid,inf);
    str = char(raw');
    meta_json = jsondecode(str);
    fclose(fid);
    
    meta_json.FS = str2double(meta_json.frameRate(1:2)); 
    
v = VideoWriter('MergedVideo.avi');% Create new video file
v.FrameRate =meta_json.FS; 
open(v)

for ii = 1:length(v_list)
    
    Vid1 = VideoReader(v_list(ii).name);
    
    while hasFrame(Vid1)
        frame = rgb2gray(readFrame(Vid1)); % read each frame
        writeVideo(v,frame) % write each frame
    end
    disp(ii)
end
close(v)