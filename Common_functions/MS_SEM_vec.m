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

% % nan_idx = isnan(data_in);
% for ii = size(data_in,2):-1:1
% 
% sem(ii) = std(data_in(:,ii), [], 1, 'omitmissing')./ sqrt(size(data_in(:,ii)));
% 
% end

sem = (std(data_in,0,1, 'omitmissing')./sqrt(size(data_in,1)));
