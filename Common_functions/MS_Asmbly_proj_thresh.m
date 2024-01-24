function [R_thresh, shuff_mat, shuff_proj] = MS_Asmbly_proj_thresh(data_in, Temp_in, nShuff, p_thresh)

if nargin < 3
    nShuff = 500;
    p_thresh = 99; 
elseif nargin < 4
    p_thresh = 99;
end



%%
rng(123,'twister')
shuff_mat = [];

for iS = 1:nShuff
    tic
    shuff_data = NaN(size(data_in));
    for ic = 1:size(data_in,2)
        
        vals = (size(data_in,1)-1).*rand(1,1) + 1;

        shuff_data(:,ic) = circshift(data_in(:,ic), floor(vals));
    end
    
    shuff_proj = assembly_activity(Temp_in,shuff_data');
    
    shuff_mat = [shuff_mat; shuff_proj];
    
end

R_thresh = prctile(shuff_mat(shuff_mat >0), p_thresh, 'all');