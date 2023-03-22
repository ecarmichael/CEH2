function ms = MS_append_timeStamps(ms, data_dir)
%% MS_append_timeStamps: grab timestamps from timeStamps.csv (assumes new Miniscope software)
%
%
%
%    Inputs: 
%    - ms : [struct] 'ms' structure from cnmfe pipeline
%
%
%
%    Outputs: 
%    - ms: [struct] 'ms' with the 'time' field added. 
%
%
%
%
% EC 2023-03-21   initial version 
%
%
%
%% initialize

if nargin < 2
    data_dir = cd;
end



%% grab timeStamps data and convert. 

TS = readtable([data_dir filesep 'timeStamps.csv']);

tvec = table2array(TS(:,2));
nan_idx = (table2array(TS(:,3)));
% correct for offsets if needed
% if tvec(1) ~= 0
%     tvec = tvec+abs(tvec(1));
% end
ms.time = tvec./1000; % convert to seconds

end