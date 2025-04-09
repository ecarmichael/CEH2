function struct_out = MS_append_mean_frame(struct_in, data_dir)
%% MS_append_mean_frame:
%
%
%
%    Inputs: 
%    - struct
%    
%    - data_dir: where the video files are stored. 
%
%
%    Outputs: 
%    - struct_out same as struct_in but with .mean_frame; 
%
%
%
%
% EC 2025-04-03   initial version 
%
%
%
%% initialize

if nargin

og_dir = cd; 


cd(data_dir)