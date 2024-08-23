function [data_out] = MS_Asmbly_participation(Proj_in, Temp_in, data_in, R_thresh)
%% MS_Asmbly_participation:
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
%  2024-07-02   initial version 
%
%
%
%% initialize

disp('bork')
data_out = NaN; 

%% loop over assemblies and get the 

for iA = size(Temp_in,2):-1:1
    
    % grab the 'significant cells in the assembly
    z_weight = zscore(Temp_in(:,iA));
    
    M_idx = z_weight > 1.96; % pull out the indicies of the best cells. 
    
    % identify significant reactivations.
    R_idx = find(Proj_in(iA,:) > R_thresh);
    
    % check the data for the cell by cell activity across each reactivation
    A_part= [];  
    
    if ~isempty(R_idx)
        for iR = length(R_idx):-1:1
            
            this_R = data_in(R_idx(iR),:)>0';
            this_R(this_R & M_idx') =2; 
            A_part(iR,:) = this_R'; 
        end
        figure(iA)
        subplot(2,1,1)
        imagesc(z_weight');
        
        subplot(2,1,2)
                    imagesc(A_part)

                    

    end
    
    
    
    
end



