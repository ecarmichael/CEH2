function [S, evts] = HF_summary_peth(S, evts)
%% HF_summary_peth: load all of the spike files and events to create summary PETHs. Meant to be quick but has limited flexibility
%
%
%
%    Inputs: 
%    - data_dir: [path] directory with the spike and event data
%
%    - cell_id: [string or cell array of strings] names of specific cells
%    to process. example HF_summary_peth(cd, {'TT2_SS_01', 'TT4_SS02})
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2025-06-10   initial version 
%
%
%
%% initialize

if nargin < 1
    error('needs both a spike ''S'' and evts ''evts'' input')
end

%%  load a

% plot check
PlotSpikeRaster([],S)
hold on
vline(evts.t{4}, 'r')


