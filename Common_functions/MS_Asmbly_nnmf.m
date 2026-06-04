function [W, H, D, opts] = MS_Asmbly_nnmf(data_in, k, opts)
%% assembly wrapper for nnmf function: 


if nargin == 1
    k = size(data_in,2);
    opts = statset('MaxIter', 10, 'Display', 'off'); % Set default options if opts is empty
elseif nargin == 2
    opts = statset('MaxIter', 10, 'Display', 'off'); % Set default options if opts is empty

end

% Initialize W and H using random values
for ii = k:-1:1
[W, H, D(ii)] = nnmf(data_in, ii, 'replicates', 10, 'Options',opts,'Algorithm','mult');
end