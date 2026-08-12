function [S, g_idx] =  MS_spike_rates(S, tvec, gap_thr)
%% appends the spike rate for each cell to the spike ts.  Based off AddRateTS from vandermeerlab codebase (OG by aacarey)


%   INPUTS
%      TS  [struct] TS structure which must contain .t with spike times
%      tvec[vector] time vector. 
%      
%      optional 
%      gap_thresh [double]  threshold for gaps in data to be excluded.  If
%      empty the data will be treated as continuous. 
%

if nargin < 3
    gap_thr = 0; 
end


% Get total time
total_time = tvec(end)-tvec(1);

% Handle gaps in time, if requested
if gap_thr > 0

    tDiff = diff(tvec);

    g_idx = find(tDiff >= gap_thr);
            
    
    LostTime = sum(tvec(g_idx + 1) - tvec(g_idx));
    total_time = total_time - LostTime;     
end

% Calculate rate and add usr field
if iscell(S.usr)
    for iS = 1:length(S.t)
        S.usr{iS}.rate = length(S.t{iS})/total_time;
    end
else
    for iS = 1:length(S.t)
        S.usr.rate(iS) = length(S.t{iS})/total_time;
    end
end

end

