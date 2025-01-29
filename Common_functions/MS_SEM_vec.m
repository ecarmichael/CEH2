function sem =  MS_SEM_vec(data_in)
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

sem = nanstd(data_in, [], 1)./ sqrt(length(data_in(~nan_idx))); 