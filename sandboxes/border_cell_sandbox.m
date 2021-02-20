%% border score sandbox


%% gen some fake cell


b_cell = zeros(10, 10);

% b_cell(1,2:8) = [2:2:6, 8,9, 4:-2:2];
b_cell(2,2:8) = [2:2:6, 8,9, 4:-2:2]*2;

% b_cell = imgaussfilt(b_cell, 1);
b_cell(b_cell < 1) = 0;

figure(20)
subplot(2,1,1)
imagesc(b_cell)


kernel_size = [size(b_cell,1) size(b_cell,2)];
 occupancy_std = 1;
[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
 Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
 kernel = pdf('Normal', Rgrid, 0, occupancy_std);
 kernel = kernel./sum(sum(kernel));
 
 b_cell_s = conv2(b_cell, kernel, 'same'); 
subplot(2,1,2)
imagesc(b_cell_s)
%% compute the border score based on solstad 2008

wall_id = {'N', 'E', 'S', 'W'}; 
   N =  sum(b_cell(1, 1:size(b_cell,2)) ~= 0)/size(b_cell,2);
   E =  sum(b_cell(1:size(b_cell,2)) ~= 0)/size(b_cell,1);
   S =  sum(b_cell(size(b_cell,2),1:size(b_cell,1)) ~= 0)/size(b_cell,1);
   W =  sum(b_cell(size(b_cell,2),1) ~= 0)/size(b_cell,2);
   
   [cM, idx] = max([N, E, S, W]); 
   
   % get the distance from the wall 
   act_bins = find((b_cell ~=0)); 
   [act_i, act_j] = find((b_cell ~=0)); 
   
   range_mat = reshape(1:(size(b_cell,1) * size(b_cell,2)), size(b_cell)); 
   
for iB = 1:length(act_bins)
    [ii, jj] = find(range_mat == act_bins(iB))
    
    
    distance=norm(cross(v1-v2,pt-v2))/norm(v1-v2)
    
end