function [start_idx, end_idx, times] = MS_get_sleep_idx(pre_post, type); 
%% MS_get_sleep_idx: extract the indicies for sleep types in concatenated data.  This is a helper function for data that has already been processed as 'all_*_post*.mat' from MS_segment_raw.m
%
%
%
%
%    Inputs: 
%    -  ms_in [struct] The 'ms_seg_resize' file from MS_segment_raw.m
%
%    - seg_idx [1 x nSegments] The 'all_Seg_idx.mat' output from
%    MS_segment_eaw.m
%
%    - type [string]  can be 'REM', 'SW' or empty for all. 
%
%
%    Outputs: 
%    - idx the start and stop indices for each segment which gets concatenated
%    together. 
%
%
%
%
% EC 2021-02-28   initial version 
%
%
%
%% initialize

if nargin < 2
    type = 'all';
end
warning off
load('ms_resize.mat', 'ms_seg_resize'); 
load('ms_trk.mat', 'ms_trk');
%% get all the start and stop indicies

   these_segs = find(contains(lower(ms_seg_resize.pre_post), lower(pre_post))); 
    

iCount = 0; 
blocks = []; 
start_idx = []; 
end_idx = []; 
for iSeg = these_segs
    if strcmpi(ms_seg_resize.hypno_label{iSeg}, type)
        iCount = iCount+1; 
        end_idx(iCount) = size(ms_seg_resize.FiltTraces{iSeg},1); 
        blocks = [blocks;  ms_seg_resize.FiltTraces{iSeg}];
        times(iCount) = ((24*60)*datenum(ms_seg_resize.time_labels{iSeg}, 'HH:MM:SS')) - ((24*60)*datenum(ms_trk.time_labels, 'HH:MM:SS')) - (ms_trk.time(end)/1000)/60;  
        fprintf([ms_seg_resize.hypno_label{iSeg} ' %i\n'], iSeg)
    end
%     disp(ms_in.pre_post{iSeg})
%     disp(seg_idx(iSeg,:))
end
for ii = length(end_idx):-1:2
    start_idx(ii) = end_idx(ii-1);
end

start_idx = start_idx(2:end) + start_idx(1:end-1);  % convert to actual indicies rather than just the length of the block

start_idx = [1 start_idx+1]; 
end_idx = [start_idx(2:end)-1 length(blocks)]; 
%% plot for sanity
% figure(10101)
% hold on
% for ii  = 1:10
%    plot(blocks(:,ii)+(ii*.1))
% end
% vline(start_idx , 'g');

