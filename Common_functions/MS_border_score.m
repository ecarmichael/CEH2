function [score, wall, wall_idx, mat_out]= MS_border_score(mat_in, smooth, smooth_sd)
%%  MS_border_score:  generates a border score for an activity matrix (mat_in) using the Solstad et al. 2008 method b = (cM -dM) / (cM + dM)
% https://science.sciencemag.org/content/322/5909/1865/tab-figures-data%
%
%    Inputs:
%    - mat_in   [n x m]  activity or significance matrix.  Ideally this is
%    a rate map or the likelihood/sig map from MI
%
%
%
%    Outputs:
%    - score  [double]  single value between -1 (center) and 1 (border
%    cell).  a score of 0.5 is typical in the lit for classifying as a border cell.  Should also check for split half stability.
%
%    - wall   [string] which wall had the lowest distance measure.  Assumes
%    N is first row of mat_in, E is last column, S is last row, and W is
%    first col.
%
%    - wall_idx  [nIndices]  indicies of the matrix corresponding to the
%    main border wall.
%
%
%
% EC 2021-02-20   initial version
%
%
%
%% initialize

if size(mat_in, 3) > 1
    error('Can only deal with 2d data ATM')
end

if nargin == 1
    smooth = 0;
    smooth_sd = 1;
elseif nargin == 2
    smooth_sd = 2; %default occupancy smoothing SD
end

%% smooth data if needed
% mat_in = abs(mat_in); % this should be done when giving a significance
% map since that contains p-values. Can be done to the data before input to
% the script. 

if smooth
%     kernel_size = [size(mat_in,1) size(mat_in,2)];
%     occupancy_std = smooth_sd;
%     [Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
%     Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
%     kernel = pdf('Normal', Rgrid, 0, occupancy_std);
%     kernel = kernel./sum(sum(kernel));
%     
%     mat_in = conv2(mat_in, kernel, 'same');
    
    mat_in = imgaussfilt(mat_in,smooth_sd);

    mat_in(mat_in < max(mat_in, [], 'all')*.3) = 0;
    
end

%% compute the border score
wall_vals = reshape(1:numel(mat_in), size(mat_in))';


wall_id = {'N', 'E', 'S', 'W'};
wall_score(1,:) =  sum(mat_in(1, 1:size(mat_in,2)) ~= 0)/size(mat_in,2);
wall_score(2,:) =  sum(mat_in(1:size(mat_in,2),size(mat_in, 1)) ~= 0)/size(mat_in,1);
wall_score(3,:) =  sum(mat_in(size(mat_in,2),1:size(mat_in,1)) ~= 0)/size(mat_in,1);
wall_score(4,:)=  sum(mat_in(1:size(mat_in,2),size(mat_in, 1)) ~= 0)/size(mat_in,2);

% get the all indicies
wall_idx(1,:) = wall_vals(1, 1:size(mat_in,2));
wall_idx(2,:) = wall_vals(1:size(mat_in,2),size(mat_in, 1));
wall_idx(3,:) = wall_vals(size(mat_in,2),1:size(mat_in,1));
wall_idx(4,:) = wall_vals(1:size(mat_in,2),size(mat_in, 1));


% get the distance from the wall
act_bins = find((mat_in ~=0)');
[act_i, act_j] = find((mat_in ~=0)');

range_mat = reshape(1:(size(mat_in,1) * size(mat_in,2)), size(mat_in))';

peak = max(mat_in,[], 'all');
[peak_idx_i, peak_idx_j] = find(mat_in == peak);

dist = [];
for iB = length(act_bins):-1:1 % look across active pixels
    for iW = size(wall_idx,1):-1:1
        if contains(wall_id(iW), {'N' 'S'})
            [ii, jj] = find(wall_vals == wall_idx(iW,:));
        else
            [ii, jj] = find(wall_vals == wall_idx(iW,:)');
        end
        [ii_pt, jj_pt] = find(range_mat == act_bins(iB));
        dist(:,iB, iW) = sqrt(((ii - ii_pt).^2)+(jj - jj_pt).^2);
    end
end

% get the distance to the closest wall
dist_m = min(dist,[], 1);

[~, best_wall] = min(sum(dist_m,2));
fprintf('Best wall is <strong>%s</strong>\n', wall_id{best_wall})

% norm to largest possible distance
dM = mean(squeeze(dist_m(:,:,best_wall)./(min(size(mat_in))/2))); % normalize to max wall length and get mean.

% get the number of active pixels on the closest wall

cM = max(wall_score);

score = (cM-dM)/(cM+dM);
fprintf('<strong>Border score = %0.2f</strong>\n', score)

wall = wall_id{best_wall};
wall_idx = best_wall;

mat_out = mat_in; % return the matrix with any smoothing that occured. 
