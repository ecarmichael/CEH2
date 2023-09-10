function sem =  MS_SEM(data_in)
%% MS_SEM:
%
%
%
%    Inputs: 
%    -
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2023-09-09   initial version 
%
%
%
%% initialize

nan_idx = isnan(data_in);

sem = std(data_in(~nan_idx)) / sqrt(length(data_in(~nan_idx))); 