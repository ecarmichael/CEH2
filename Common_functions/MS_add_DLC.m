function MS_add_DLC(behav_in, DLC)
%% MS_add_DLC: appends DeepLabCut position data to the ms behav.mat struct for each measure in the DLC.
%
%
%    Inputs: 
%    -  behav_in: [struct]  contains time, positions,...
%
%    -  DLC:  [ntime x nParts] OR [string]  output from DLC as a .csv file.
%    A file name can be given as an input which will load and process the
%    raw file. 
%
%    Outputs: 
%    -
%
%
%
%
% EC 2021-02-17   initial version 
%
%
%
%% initialize