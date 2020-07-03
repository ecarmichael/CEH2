function ms = MS_msExtractBinary_detrendTraces(ms, z_threshold)
%MSEXTRACTBINARY Summary of this function goes here
%   Detailed explanation goes here
%
%
%  EC update June 2020:   added threshold as an argin. 

%% check inputs

if nargin ==1
    z_threshold = 3;
    fprintf('<strong>%s: </strong> No threshold specified by user. Used Z-score of 3\n', mfilename);
else
    fprintf('<strong>%s: </strong> User specified Z-score = %d\n', mfilename, z_threshold);
end

%% Parameters

    [bFilt,aFilt] = butter(2,  2/(30/2), 'low');

    for trace_i = 1:ms.numNeurons
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

