function ms = MS_msExtractBinary_detrendTraces(ms, z_threshold)
% MS_msExtractBinary_detrendTraces  EC version of
% msExtractBinary_detrendTraces with added zthreshold and made the sampling
% frequency in the butterworth filter design adaptive to ms.time using
% mode(diff(ms.time)). 
%
%
%  EC update June 2020:   added threshold as an argin. 
%  EC update July 2020:   added adaptive Fs in butter. 

%% check inputs

if nargin ==1
    z_threshold = 3;
    fprintf('<strong>%s: </strong> No threshold specified by user. Used Z-score of 3\n', mfilename);

else
    fprintf('<strong>%s: </strong> User specified Z-score = %d\n', mfilename, z_threshold);
end

%% Parameters
    fprintf('<strong>%s: </strong> Infered sampling rate = %d Hz\n', mfilename, floor(1/(0.001*(mode(diff((ms.time)))))));

    [bFilt,aFilt] = butter(2,  2/(floor(1/(mode(diff(ms.time))))/2), 'low');

    for trace_i = 1:size(ms.FiltTraces,2)
        detrend_raw_trace=detrend(ms.FiltTraces(:,trace_i),2);   
        filt_trace = zscore(filtfilt(bFilt,aFilt,detrend_raw_trace));

        diff_trace = diff(filt_trace);
        diff_trace(end+1) = 0;

        binary_trace = filt_trace*0;
        binary_trace(filt_trace>z_threshold & diff_trace>0) = 1;

        ms.detrendRaw(:,trace_i)=detrend_raw_trace;
        ms.Binary(:,trace_i) = binary_trace;

    end

ms.Binary_threshold = z_threshold;

end

