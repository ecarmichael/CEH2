function h = MS_gscatter3d(x,y,z, idx, c)
%% MS_gscatter3d: wrapper to make the built-in gscatter function take in z data. 
%
%   - MS_gscatter3d(x,y,z,...): x, y, z are data vectors of the same
%   length.
%
%   - MS_gscatter3d(X, ...): X is a [3 x n] array. 
%
%
%    Inputs: 
%     - x, y, z are [1 x n] vetors of the same length.  
%     - X [3 x n] array
%     - idx: [1 x n] vector of grouping labels (ex: output from kmeans).
%     - c: [3 x n] color code.  [ONLY works with values not STRINGS]
%
%
%    Outputs: 
%     - h : figure handle
%
%
%
%
% EC 2020-02-17   initial version based on 'gscatter' Copyright 1993-2019 The MathWorks, Inc.
%
%
%%  inititalize

if size(x,1) ==3
    x = x';
end

% too few
if nargin < 2
    error('<strong>MS_gscatter3d</strong>: Too few inputs')
    
    % using X [3 x N] but no c specified. 
elseif nargin ==2 && size(x,2) ==3
    fprintf('<strong>MS_gscatter3d</strong>: using X [%d x %d] but no c specified.\n', size(x,2), size(x,1))
    method = 'X';
    idx = y; % swap spots. 
    clear y
    c = linspecer(3); 
    
elseif nargin ==3 && size(x,2) ==3
    fprintf('<strong>MS_gscatter3d</strong>: using X [%d x %d].\n', size(x,2), size(x,1))
    method = 'X';
    idx = y;
    clear y; 
    c = z; 
    
elseif nargin ==3 && size(x,1) ~=3
    error('<strong>MS_gscatter3d</strong>: invalid input for x. Must be either 3 x n or n x 3, or 1 x n')
    
elseif nargin <= 3 && nargin ~=5
    error('<strong>MS_gscatter3d</strong>: inputs must be either MS_gscatter3d(x,y,z,idx,c) or MS_gscatter3d(X,idx, c)')
    
elseif nargin ==4 && size(x,1) == size(y,1) && size(x,1) == size(z,1) && size(y,1) == size(z,1)
    fprintf('<strong>MS_gscatter3d</strong>: X input is %d x %d. Using x, y, z as 1 x n vectors of equal length. no c specified.\n', size(x,2), size(x,1))
    method = 'xyz';
    c = linspecer(3); 
    
elseif nargin ==5 && size(x,1) == size(y,1) && size(x,1) == size(z,1) && size(y,1) == size(z,1)
    fprintf('<strong>MS_gscatter3d</strong>: X input is %d x %d. Using x, y, z as 1 x n vectors of equal length.\n', size(x,2), size(x,1))
    method = 'xyz';
end
    

%% plot
switch method
    
    case 'xyz'

        h = gscatter(x,y,idx,c);
        hold on
        idx_u = unique(idx);
        for ii = 1:numel(idx_u)
            set(h(ii), 'ZData', z(idx == idx_u(ii),1));
        end

    case 'X'
        
        h = gscatter(x(:,1), x(:,2),idx,c);
        hold on
        idx_u = unique(idx);
        for ii = 1:numel(idx_u)
            set(h(ii), 'ZData', x(idx == idx_u(ii),3));
        end

end

view(3)
