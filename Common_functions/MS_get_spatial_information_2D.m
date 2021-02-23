function [MI, posterior, occupancy_vector, p_active, likelihood] = MS_get_spatial_information_2D(binary_in, position_in, X_bin_vec, Y_bin_vec)
%% MS_get_spatial_information: extracts the spatial information for a given binary signal across positions. 
%
%  This script is based on the one by Guillaume Etter (2017-) 
%  @ https://github.com/etterguillaume/CaImDecoding/blob/master/extract_2D_information.m
%  and can be credited to this publication: 
%  Etter G, Manseau F and Williams S (2020) 
%  A Probabilistic Framework for Decoding Behavior From in vivo Calcium Imaging Data. 
%  Front. Neural Circuits 14:19. doi: 10.3389/fncir.2020.00019
%
%  This version is only a modification to fit the MS codebase. It was an
%  exercise in understanding the methods used in the original script. 
%
%    Inputs: 
%    - binary_in [1 x nSamples] binary vector from a Ca trace.  Binary
%    thresholds are set using MS_msExtractBinary_detrendTraces
%
%    - binary_tvec [1 x nSamples] time vector corresponding to the Ca trace
%
%    - position_in [1 or 2 x nSamples] position array for X and Y over
%    time. 
%
%    Outputs: 
%    - MI:   mutual information between the position and binary activity
%
%    - posterior:  
%
%
%
% EC 2020-07-17   initial version 
%
%
%
%% initialize

if nargin <4
    error('requires binary_in, position_in(n x2), X bins, Y bins, as inputs')
end




%% Create bin vectors
p_active = sum(binary_in)./length(binary_in); % Expressed in probability of firing (<1)

%% Compute joint probabilities (of cell being active while being in a state bin)
likelihood = zeros(length(X_bin_vec)-1,length(Y_bin_vec)-1);
occupancy_vector = zeros(length(X_bin_vec)-1,length(Y_bin_vec)-1);
MI = 0;

for iY = 1:length(Y_bin_vec)-1
    for iX = 1:length(X_bin_vec)-1
    binarized_spatial_vector = 0*binary_in;
    position_idx = find(position_in(:,1) >= X_bin_vec(iX) & position_in(:,1) < X_bin_vec(iX+1) & position_in(:,2) >= Y_bin_vec(iY) & position_in(:,2) < Y_bin_vec(iY+1));
    
    if ~isempty(position_idx)
        binarized_spatial_vector(position_idx)=1;
        occupancy_vector(iX, iY) = length(position_idx)/length(binary_in);
        activity_in_bin_idx = find(binary_in == 1 & binarized_spatial_vector == 1);
        inactivity_in_bin_idx = find(binary_in == 0 & binarized_spatial_vector == 1);
        likelihood(iX, iY) = length(activity_in_bin_idx)/length(position_idx);
        
        joint_prob_active = length(activity_in_bin_idx)./length(binary_in);
        joint_prob_inactive = length(inactivity_in_bin_idx)./length(binary_in);
        prob_in_bin = length(position_idx)./length(binary_in);
        
        if joint_prob_active ~= 0
            MI = MI + joint_prob_active*log2(joint_prob_active./(prob_in_bin*p_active));
        end
        if joint_prob_inactive ~= 0
            MI = MI + joint_prob_inactive*log2(joint_prob_inactive./(prob_in_bin*(1-p_active)));
        end
    end
    end
end

posterior = likelihood.*occupancy_vector/p_active;


