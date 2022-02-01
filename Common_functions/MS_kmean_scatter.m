function [g_idx, n_val, n_freq] = MS_kmean_scatter(data_in, nGroups, dim, markersize)
%% kmeans 3d scatter with grouped colours 



% EC 2020 initial version based on https://www.mathworks.com/matlabcentral/answers/225193-grouping-data-for-a-3d-scatter-plot


%% initialize

if nargin <2
    error('Not enough inputs: requires data_in and nGroups')
elseif nargin < 3
    dim = 1:length(data_in);
    markersize = 20; 
elseif nargin < 4
    markersize = 20; 

end


%% compute kmeans

[idx, ~] = kmeans(data_in, nGroups); 
u=unique(idx);                % find the unique values in vectr
[n,b]=histc(idx,u);           % bin on those and count number in bin
[n_val,is]=sort(n,'descend');   % sort highest first and keep index array of position
 n_freq=idx(arrayfun(@(x) find(b==x,1,'first'),is));

%% plot
[uni_groups, ~, IC] = unique(idx); % get the index values 
colors = linspecer(length(uni_groups));  %or any other way of creating the colormap
% markersize = 20;   %change to suit taste
scatter3(data_in(:,dim(1)), data_in(:,dim(2)), data_in(:,dim(3)), markersize, colors(IC,:), 'filled');

%% output the grouping index

g_idx = idx; 
