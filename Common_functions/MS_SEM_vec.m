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
for ii = length(data_in):-1:1

sem(ii) = nanstd(data_in(:,ii), [], 1)./ sqrt(sum(~isnan(data_in(:,ii))));

end