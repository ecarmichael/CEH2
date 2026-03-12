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

% nan_idx = isnan(data_in);
for ii = size(data_in,2):-1:1

sem(ii) = std(data_in(:,ii), [], 1, 'omitmissing')./ sqrt(sum(data_in(:,ii), 'omitmissing'));

end