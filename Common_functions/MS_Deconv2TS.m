function TS = MS_Deconv2TS(cfg_in, ms)
%% MS_Binary2TS: converts a binary Ca trace into timestamps 'TS'
%
%
%
%    Inputs: 
%     - cfg: [struct] configuration (see defaults below. 
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
%% set up defaults
% cfg_def.units = 'sec'; % convert default ms.time (ms) to seconds. 

% cfg = ProcessConfig(cfg_def, cfg_in); 
%% convert

 Csp = ms.deconv./ms.denoise; 
 Csp = Csp > 0.01;

for iC = size(ms.deconv, 2):-1:1
    timestamps{iC} = ms.time(diff(Csp(:,iC)) > 0);
    label{iC} = num2str(iC); 
%     switch cfg.units
%     case 'sec'
%         timestamps{iC} = timestamps{iC}/1000;
%     case 'ms'
%         timestamps{iC} = timestamps{iC};
%     end
%     centroid{iC} = ms.Cen  [todo: figure out what is going on with the
%     inconsistent numbers here. 
end


TS = ts(timestamps);
TS.label = label;
TS.units = 'sec';

