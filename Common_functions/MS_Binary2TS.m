function TS = MS_Binary2TS(ms)
%% MS_Binary2TS: converts a binary Ca trace into timestamps 'TS'
%
%
%
%    Inputs: 
%     - ms:  [struct] common miniscope ms struct.  Requires Binary field from msExtractBinary_detrendTraces
%
%
%
%    Outputs: 
%     - TS  [struct]  'TS' data struct containing the timestamps for each
%     binary transition.  
%
%
%
%
% EC 2020-04-24   initial version 
%
%
%% convert
for iC = size(ms.Binary, 2):-1:1
    timestamps{iC} = ms.time(diff(ms.Binary(:,iC)) > 0);
    label{iC} = num2str(iC); 
end

TS = ts(timestamps);
TS.label = label;

