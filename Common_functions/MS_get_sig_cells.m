function [p_sig, sig_vec] = MS_get_sig_cells(All_cells, threshold, metric)
%% MS_get_sig_cells:  finds all cells that meet the threshold value
%
%
%
%    Inputs: 
%    - All_cells  [struct] contains the output from Spatial_screener_info. 
%
%    - field [string] which field to process.  Eg 'Place', 'Speed', 'Acc'
%
%    - threshold:  pvalue threshold
%
%    Outputs: 
%    - p_sig:  percentage of significant 
%
%    - sig_vec:  [1 x nCell]  Logical vector os significance for each cell.
%    
%
%
% EC 2021-04-27   initial version 
%
%
%
%% initialize


if nargin == 0
    load('SI.mat', 'SI');
    threshold = 0.01;  
elseif nargin == 1
    threshold = 0.01;

end

% fprintf('<strong> %s </strong>: p-value threshold set at <strong>%0.2f</strong>\n', mfilename, threshold)

%% cycle through cells and check sig for this field. 
if ~exist('metric','var')
    fnames = fieldnames(All_cells(1).spatial);
else
    fnames{1} = metric;
end

sig_vec = zeros(length(All_cells),length(fnames));
for iF = 1:length(fnames)
    if ~strcmp(fnames{iF}, 'fname') && ~strcmp(fnames{iF}, 'events')

        for iC = length(All_cells):-1:1
            if (All_cells(iC).spatial.(fnames{iF}).MI_pval < threshold) && ((All_cells(iC).spatial.(fnames{iF}).MI_pval) ~= 0)
                sig_vec(iC,iF) = 1;
            end
        end
        p_sig(iF) = sum(sig_vec(:,iF))/length(sig_vec(:,iF));
        fprintf('<strong>%0.0f/%0.0f</strong> cells (%0.2f%%) passed significance threshold of %0.2f for <strong>%s</strong>\n',sum(sig_vec(:,iF)), length(sig_vec(:,iF)),p_sig(iF)*100, threshold, fnames{iF})
    end
end





