function [ReAct_stats, shuff_mat, shuff_proj, all_shuff_proj] = MS_Asmbly_proj_thresh(data_in, Temp_in, nShuff, p_thresh, Proj_in)

if nargin < 3
    nShuff = 1000;
    p_thresh = 99;
    Proj_in = [];
elseif nargin < 4
    p_thresh = 99;
    Proj_in = [];
elseif nargin < 5
    Proj_in = [];
end

if isempty(Temp_in)
    ReAct_stats.R_thresh = NaN;
    shuff_proj = [];
    shuff_mat = [];
    return
end

%%
rng(123,'twister')
shuff_mat = [];
all_shuff_proj = []; 
for iS = nShuff:-1:1
    tic
    shuff_data = NaN(size(data_in));
    for ic = 1:size(data_in,2)
        
        c_idx = randsample(1:size(data_in,1), 1);

        shuff_data(:,ic) = circshift(data_in(:,ic), c_idx);
    end
    
    shuff_proj = assembly_activity(Temp_in,shuff_data');
    % hold all the projections as a way to check projection stats. 
    all_shuff_proj{iS} = shuff_proj; 

    try
        this_idx = randsample(1:size(shuff_proj,1), 1);
    catch 
        display(size(shuff_proj))
    end

        shuff_mat(iS,:) = shuff_proj(this_idx,:);

end

ReAct_stats.R_thresh = prctile(shuff_mat(shuff_mat >0), p_thresh, 'all');

end