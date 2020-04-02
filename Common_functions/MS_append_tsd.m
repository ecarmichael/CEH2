function TSD_out = MS_append_tsd(TSD_in, nargin)
%% MS_append_tsd: appends TSD structs (typically CSCs) together in the order that the exist in the TSD_in input. 
%
%  [WARNING] this is a basic function that only takes in 
%
%
%
%    Inputs: 
%     - TSD_in: [struct] a single cell array of TSDs
%     structs.  
%           ex: TSD_in{1} + TSD_in{2} + ... + TSD_in{end}  ==> TSD_out
%
%
%
%    Outputs: 
%     - a single new TSD_out
%        TSD_out:
%                  type: 'tsd'
%                 units: 'V'
%                  tvec: [nSamples×1 double]
%                  data: [nChannel×nSamples double]
%                 label: {'CSC1'  'CSC2' ...}
%                   cfg: [1×1 struct]
%
%
%
%
% EC 2020-04-01   initial version 
%
%
%% determine the input data type
if ~iscell(TSD_in)
    error(sprintf('<strong>%s</strong>: TSD_in is not a cell', mfilename));
end

%% cycle through all the TSD_in cells and append them

TSD_out = TSD_in{1};
TSD_out.data = [];
TSD_out.tvec = [];

for iT = 1:length(TSD_in)
    TSD_out.data = [TSD_out.data, TSD_in{iT}.data];
    TSD_out.tvec = [TSD_out.tvec; TSD_in{iT}.tvec];

end

