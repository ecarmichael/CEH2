function [P_temp, P_proj, P_cells] = MS_Asmbly_select(A_temp, A_proj, threshold)
%% MS_Asmbly_select:
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
% EC 2024-01-24   initial version 
%
%
%
%% initialize

keep_idx = zeros(1,size(A_temp,2));
P_cells = cell(size(keep_idx));
for ii = size(A_temp,2):-1:1
    
    if max(A_temp(:, ii)) > 0.2
        keep_idx(ii) = 1;
        
        z_weight = zscore(A_temp(:,ii));
        
        P_cells{ii} = find(z_weight > 1.96);
        
    end

end

P_temp = A_temp;
P_temp(:,~keep_idx) = [];

P_cells(~keep_idx) = []; 

P_proj = A_proj;
P_proj(~keep_idx,:) = [];