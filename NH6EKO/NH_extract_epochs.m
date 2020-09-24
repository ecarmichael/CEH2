function [start_idx, end_idx, tran_states]= NH_extract_epochs(data)
%% NH_extract_epochs: scans through the data variable to identify transitions from one value to another and extract the length of blocks.
%
%
%
%
%    Inputs:
%    - data:  [1 x nSamples] vector contaning some values to be broken down into blocks.
%
%    [optional]
%    - tvec:  [1 x nSamples]  time vector which can be used for computing
%    durations. must be the same length as data
%
%    Outputs:
%    - start_idx:  {nUniqueValues} continaing [1 x nBlocks]  start index of each block (cells)
%    - end_idx:  {nUniqueValues} continaing [1 x nBlocks]  end index of each block (cells)
%    - trans_states:  {nUniqueValues} continaing [1 x nBlocks]  the value of
%    the next block as a number where the first value is the current state
%    and the next is the subsequent state.  example: trans_states(1) =  32
%    would be a transition from state '3' to state '2'. 
%
%     Example: data = [ 1 1 1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 1 1 1 1 3 3 3 ];
%         start_idx{1} = [1 17]
%         end_idx{1} = [7 24]
%         trans_states{1} = [12, 13]
% EC 2020-09-17   initial version
%
%
%
%% initialize
if nargin ==0
    error('Requires at least one input')
end
%% determine the values in the score input
tic; uVals = unique(data); toc; % this is slow

uVals(isnan(uVals)) = []; % remove any NaNs which for some reason get multiple counts in unique
[mrows, ncols] = size(uVals');
outputstr = ['%' num2str(mrows) 'i ']; % template for the string, you put your datatype here
outputstr = repmat(outputstr, 1, ncols); % replicate it to match the number of columns
outputstr = [outputstr '\n']; % add a new line if you want
fprintf(['<strong>%s</strong>: found %i unique values: <strong>' outputstr '</strong> \n'], mfilename, length(uVals), uVals)

%% cycle through each value and extract both the nSamples per block and the transition

for iVal = length(uVals):-1:1
    diff_val = diff(data == uVals(iVal));
    start_idx{iVal} = find(diff_val == 1);
    end_idx{iVal} = find(diff_val == -1);
    
    % corrections
    if start_idx{iVal}(end) > end_idx{iVal}(end) % correct for cases where it ends with the desired value.
        end_idx{iVal} = [end_idx{iVal}; length(data)];
        tran_val{iVal} = data(end_idx{iVal}(1:end-1)); %transition values
    else
        tran_val{iVal} = data(end_idx{iVal}+1); %transition values
    end
    if start_idx{iVal}(1) > end_idx{iVal}(1) % correct for cases where it starts with the desired value.
        start_idx{iVal} = [1; start_idx{iVal}];
    end
    for ii = length(tran_val{iVal}):-1:1
        tran_states{iVal}(ii) = str2double([num2str(iVal) num2str(tran_val{iVal}(ii))]); 
    end
    
end




