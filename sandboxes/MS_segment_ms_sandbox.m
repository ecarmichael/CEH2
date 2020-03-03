function ms_out = MS_segment_ms_sandbox(cfg_in,ms_in)
%% MS_segment_ms: segments the previously concatinated data in the ms stucture
%   MS_segment_ms takes the previously concatenated and processed miniscope
%   data in the 'ms' data structure from
%   'msGenerateVideoObj_ConcatenateBatch' in
%   https://github.com/torilynntemple/MiniscopeAnalysisConcatenation by
%   Tori-Lynn Temple (Williams lab at McGill) and breaks it back down into
%   segments based on unique recording periods.  
%
%   Inputs:
%       ms_in: [struct] this is the data structure as output by msGenerateVideoObj_ConcatenateBatch
%       
%   Outputs:
%       ms_out: [struct] same as ms_in but with new segmented fields for
%       time, RawTraces, FiltTraces, vidNum, frameNum. 
%
%
%   To Do:
%       - get file name, date, and other identifiers for the video files
%       - have a feature that checks the camera idea based on the video
%       file names (behav vs mscope) and check the samples add up to frames
%       - append LFP data (provided alignment is spot on...
%       - figure out what is happeneing with the ms.timestamps0 and
%       ms.timestamps1 being different. missing frame correction? 
%
%
%
% EC: 2019-12-13:
%   - initial sandbox with basic segmentation using samples in the
%   ms.timestamps and ms.time fields
%% intialize
cfg_def = [];
cfg_def.user_fields = {}; % user fields to remove: example ' BinaryTraces'

cfg = ProcessConfig(cfg_def, cfg_in); 



ms_out = ms_in; 

% remove the fields that will be replaced (for cleanliness atm)
ms_out = rmfield(ms_out, 'RawTraces');
ms_out = rmfield(ms_out, 'FiltTraces');
ms_out = rmfield(ms_out, 'time');
ms_out = rmfield(ms_out, 'vidNum');
ms_out = rmfield(ms_out, 'frameNum');

%user specified fields. 
if ~isempty(cfg.user_fields)
    for iF = 1:length(cfg.user_fields)
        ms_out = rmfield(ms_out, cfg.user_fields{iF});
    end
end


%% get all the time blocks and restrict the data to those blocks
times = ms_in.timestamps1;

% make incremental 
all_times = nan(1,length(times));
for iT = 1:length(times)
    if iT == 1
        all_times(iT)  =  times(iT);
    else
        all_times(iT)  =  times(iT) + all_times(iT-1);
    end
end

% check 
if all_times(end) ~= length(ms_in.time)
    error('Time segments and the total number of frames from ms.time do not align')
end
    
% loop though the times and restrict to the individual segments    
for iT = 1:length(all_times)
    if iT == 1 
        ms_out.time{iT,1} = ms_in.time(1:all_times(iT));
        ms_out.RawTraces{iT,1} = ms_in.RawTraces(1:all_times(iT),:);
        ms_out.FiltTraces{iT,1} = ms_in.FiltTraces(1:all_times(iT),:);
        if ~isempty(cfg.user_fields)
            for iF = 1:length(cfg.user_fields)
                ms_out.(cfg.user_fields{iF}){iT,1} = ms_in.(cfg.user_fields{iF})(1:all_times(iT),:);
            end
        end
        ms_out.frameNum{1,iT} = ms_in.frameNum(1,1:all_times(iT));
        ms_out.vidNum{1,iT} = ms_in.vidNum(1,1:all_times(iT));
    else 
        ms_out.time{iT,1} = ms_in.time(all_times(iT-1)+1:all_times(iT));
        ms_out.RawTraces{iT,1} = ms_in.RawTraces(all_times(iT-1)+1:all_times(iT),:);
        ms_out.FiltTraces{iT,1} = ms_in.FiltTraces(all_times(iT-1)+1:all_times(iT),:);
        if ~isempty(cfg.user_fields)
            for iF = 1:length(cfg.user_fields)
                ms_out.(cfg.user_fields{iF}){iT,1} = ms_in.(cfg.user_fields{iF})(all_times(iT-1)+1:all_times(iT),:);
            end
        end
        ms_out.frameNum{1,iT} = ms_in.frameNum(1,all_times(iT-1)+1:all_times(iT));
        ms_out.vidNum{1,iT} = ms_in.vidNum(1,all_times(iT-1)+1:all_times(iT));
    end
    if length(ms_out.time{iT,1}) ~= times(iT)
        error('Time segments are not aligning propertly.')
    end
end

%% clean up

ms_out.format = 'Segmented based on timestamp periods only'; 

disp('Segmented based on timestamp periods only')

