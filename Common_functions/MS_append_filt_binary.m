function [ms_out, Binary_out] = MS_append_filt_binary(ms, z_threshold)
%% : MS_append_filt_binary: extracts the binary 'on' phases from the denoised signal instead of the filtered signal as per msExtractBinary_detrendTraces.m
%
%
%
%    Inputs:
%    - ms: [struct]      ms structure containing the denoised signal
%
%    - z_threshold: [double]     cut off value (in std) for the zscored
%    data.
%
%    Outputs:
%    - ms: [struct]     ms structutre with the binary appended to it.
%
%
%
%
% EC 2023-06-14   initial version
%
%
%
%% initialize

ms_out = ms;

%% Parameters
if nargin < 2
    z_threshold = 2;
end

if iscell(ms.FiltTraces)
    for ii = length(ms.FiltTraces):-1:1
        
        for trace_i = ms.numNeurons:-1:1
            this_trace = ms.FiltTraces{ii}(:, trace_i);
            this_trace = zscore(this_trace);
            
            d1_trace = [diff(this_trace); 0];
            
            binary_trace = this_trace*0;
            binary_trace(this_trace>z_threshold & d1_trace>0) = 1;
            
            ms_out.Binary{ii}(:,trace_i) = binary_trace;
            
        end
    end
    
    
else
    for trace_i = ms.numNeurons:-1:1
        this_trace = ms.denoise(:, trace_i);
        this_trace = zscore(this_trace);
        
        d1_trace = [diff(this_trace); 0];
        
        
        binary_trace = this_trace*0;
        binary_trace(this_trace>z_threshold & d1_trace>0) = 1;
        
        ms_out.Binary(:,trace_i) = binary_trace;
        
    end
end

    ms_out.Binary_cfg.z_thresh = z_threshold;
    ms_out.Binary_cfg.method = 'Filt';
    
    Binary_out  = ms_out.Binary; 


