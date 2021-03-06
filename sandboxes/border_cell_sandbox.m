%% border score sandbox


%% gen some fake cell


b_cell = zeros(10, 10);

% b_cell(1,2:8) = [2:2:6, 8,9, 4:-2:2];
b_cell(1,2:8) = [2:2:6, 8,9, 4:-2:2]*4;
b_cell = b_cell./max(b_cell,[],'all'); % normalize activity. 
% b_cell = imgaussfilt(b_cell, 1);

figure(20)
subplot(2,2,1)
imagesc(b_cell)


kernel_size = [size(b_cell,1) size(b_cell,2)];
 occupancy_std = 1;
[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
 Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
 kernel = pdf('Normal', Rgrid, 0, occupancy_std);
 kernel = kernel./sum(sum(kernel));
  
 b_cell_s = conv2(b_cell, kernel, 'same');
 
  b_cell_s(b_cell_s < max(b_cell_s, [], 'all')*.3) = 0;

subplot(2,2,3)
imagesc(b_cell_s)

% generate a center cell for comparison
c_cell = zeros(10, 10);

c_cell(5:7,4:6) = [4,6,4; 6, 8, 6; 4, 6, 4];
c_cell = c_cell./max(c_cell,[],'all'); % normalize activity. 


kernel_size = [size(c_cell,1) size(c_cell,2)];
 occupancy_std = 1;
[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
 Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
 kernel = pdf('Normal', Rgrid, 0, occupancy_std);
 kernel = kernel./sum(sum(kernel));
  
 c_cell_s = conv2(c_cell, kernel, 'same');
 
 c_cell_s(c_cell_s < max(c_cell_s, [], 'all')*.3) = 0;
 
 
 % generate a corner edge cell for comparison
e_cell = zeros(10, 10);

e_cell(1:3,1:3) = [4,6,4; 6, 8, 6; 4, 6, 4];
e_cell = e_cell./max(e_cell,[],'all'); % normalize activity. 


kernel_size = [size(e_cell,1) size(e_cell,2)];
 occupancy_std = 1;
[Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
 Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
 kernel = pdf('Normal', Rgrid, 0, occupancy_std);
 kernel = kernel./sum(sum(kernel));
  
 e_cell_s = conv2(e_cell, kernel, 'same');
 
 e_cell_s(e_cell_s < max(e_cell_s, [], 'all')*.3) = 0;


%% compute the border score based on solstad 2008

% for testing
wall_vals = reshape(1:numel(b_cell_s), size(b_cell_s))';


wall_id = {'N', 'E', 'S', 'W'}; 
   wall_score(1,:) =  sum(b_cell_s(1, 1:size(b_cell_s,2)) ~= 0)/size(b_cell_s,2);
   wall_score(2,:) =  sum(b_cell_s(1:size(b_cell_s,2),size(b_cell_s, 1)) ~= 0)/size(b_cell_s,1);
   wall_score(3,:) =  sum(b_cell_s(size(b_cell_s,2),1:size(b_cell_s,1)) ~= 0)/size(b_cell_s,1);
   wall_score(4,:)=  sum(b_cell_s(1:size(b_cell_s,2),size(b_cell_s, 1)) ~= 0)/size(b_cell_s,2);
   
   % get the all indicies
   wall_idx(1,:) = wall_vals(1, 1:size(b_cell_s,2));
   wall_idx(2,:) = wall_vals(1:size(b_cell_s,2),size(b_cell_s, 1));
   wall_idx(3,:) = wall_vals(size(b_cell_s,2),1:size(b_cell_s,1));
   wall_idx(4,:) = wall_vals(1:size(b_cell_s,2),size(b_cell_s, 1));

   [cM, idx] = max(wall_score); 
   
   % get the distance from the wall 
   act_bins = find((b_cell ~=0)'); 
   [act_i, act_j] = find((b_cell ~=0)'); 
   
   range_mat = reshape(1:(size(b_cell,1) * size(b_cell,2)), size(b_cell))'; 
   
   peak = max(b_cell,[], 'all'); 
   [peak_idx_i, peak_idx_j] = find(b_cell == peak); 
   
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
 dM = mean(squeeze(dist_m(:,:,best_wall)./(min(size(b_cell_s))/2))); % normalize to max wall length and get mean. 
 
 
 b_score = (cM-dM)/(cM+dM); 
      fprintf('Border score = %d</strong>\n', b_score)

      
      %% check plot
      
figure(20)
subplot(2,3,1)
imagesc(b_cell)
b_score = MS_border_score(b_cell);
text(0, 0, ['Border score = ' num2str(b_score,2)])
    
subplot(2,3,2)
imagesc(c_cell)
b_score = MS_border_score(c_cell);
text(0, 0, ['Border score = ' num2str(b_score,2)])

subplot(2,3,3)
imagesc(e_cell)
b_score = MS_border_score(e_cell);
text(0, 0, ['Border score = ' num2str(b_score,2)])

subplot(2,3,4)
imagesc(b_cell_s)
b_score = MS_border_score(b_cell_s);
text(0, 0, ['Border score = ' num2str(b_score,2)])

subplot(2,3,5)
imagesc(c_cell_s)
b_score = MS_border_score(c_cell_s);
text(0, 0, ['Border score = ' num2str(b_score,2)])

subplot(2,3,6)
imagesc(e_cell_s)
b_score = MS_border_score(e_cell_s);
text(0, 0, ['Border score = ' num2str(b_score,2)])

%% Get borderscores for all cells

for iC = 1:length(All_cells.place.MI)
b_score = MS_border_score(abs(All_cells.place.Sig_map(:,:,iC)), 1, 1); 
All_cells.border.score(iC) = b_score;
end


figure(303)
bins = -1:.1:1; 
h = histogram(All_cells.border.score, bins) ;
h.FaceColor = PARAMS.red;  h.FaceAlpha = .9;
h.EdgeColor = 'k';  h.EdgeAlpha = 1; 
h.LineWidth = 1;
xlabel('border score')
ylabel('count')
title('border score for ck2cre-1359 2021-01-30')
vline(0.5, 'r', 'thresh')

for iC = 1:round(length(All_cells.border.score)/3)
    if All_cells.border.score(iC) >= 0.4  && sum(All_cells.place.Sig_map(:,:,iC) ~= 0, 'all') > 3
        [~,~,~,this_mat] = MS_border_score(All_cells.place.Sig_map(:,:,iC), 1, 1); 

        figure(1000+ iC)
        imagesc(0:3:30, 0:3:30, this_mat)
        axis xy
        text(0,33, ['Cell ID: ' num2str(iC) '  Border Score: ' num2str(All_cells.border.score(iC), 2)])
    end
end
    
