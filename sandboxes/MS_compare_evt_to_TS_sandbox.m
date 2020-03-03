function evt_out = MS_compare_evt_to_TS_sandbox(cfg_in,evts_in, TS_in)
%% MS_compare_evt_to_TS: used for aligning the Neuralynx .nev events to the Miniscope timestamp.dat 
%   files.  
%
%   WIP: needs ground truth data to validate. 
%
%    Inputs: 
%     - cfg_in: [struct] user configuration paramters. if empty [] then use
%     defaults (see below)
%     - evts_in: [struct] Neuralynx events in the vanderMeer lab format
%     from loadEvents
%     - TS_in: [struct] 
%
%    Outputs: 
%     -
%
%
%
%
% EC 2020-01-14   initial version 
%
%
%% cycle through the events blocks found in the NLX .nev and compare them to the ones in the TS 
% structure containing the timestamp.dat data.




disp('Compare')

for iRec = 1:length(evts_in.tstart)
    %     disp(['Rec ' num2str(iRec)])
        fprintf('Evts id: %.0f = %.0f samples fs ~ %.1f time: %0.2f sc || TS id: %.0f = %.0f samples fs ~ %.1f time: %0.2f sc\n',...
            iRec,...
            length(rec_evt{iRec}.t{this_evt}),...
            mode(diff(rec_evt{iRec}.t{this_evt}))*1000,...
            length(rec_evt{iRec}.t{this_evt})/(1/(median(diff(rec_evt{iRec}.t{this_evt})))),...
            iRec, length(TS{iRec}.system_clock{end}),...
            mode(diff(TS{iRec}.system_clock{end})),...
            length(TS{iRec}.system_clock{1})/(1/(median(diff(TS{iRec}.system_clock{1}(2:end)))*0.001)))
    end
    evt_TS_diff(iRec) = length(rec_evt{iRec}.t{this_evt}) - length(TS{iRec}.system_clock{end});
end

evt_TS_diff % print the offset






