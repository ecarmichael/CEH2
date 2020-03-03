function [ms] = msExtractBinary_detrendTraces(ms);
%MSEXTRACTBINARY Summary of this function goes here
%   Detailed explanation goes here

%% Parameters
    z_threshold = 3;

    [bFilt,aFilt] = butter(2,  2/(30/2), 'low');

    for trace_i = 1:ms.numNeurons;
        detrend_raw_trace=detrend(ms.RawTraces(:,trace_i),2);   
        filt_trace = zscore(filtfilt(bFilt,aFilt,detrend_raw_trace));

        d1_trace = diff(filt_trace);
        d1_trace(end+1) = 0;

        binary_trace = filt_trace*0;
        binary_trace(filt_trace>z_threshold & d1_trace>0) = 1;

        ms.detrendRaw(:,trace_i)=detrend_raw_trace;
        ms.Binary(:,trace_i) = binary_trace;

    end



end

