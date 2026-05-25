function [A_Maps, s_idx] = MS_Asmbly_ratemap(pos_vec, proj, bins)
%% MS_assmbly_ratemap: makes a rate map for each projection for the positions within the bins. 


%% initialize the bins

A_Maps = zeros(size(proj,1),size(bins,2));

% Calculate the rate for each bin
for ii = 1:length(bins)-1
    binMask = (pos_vec >= bins(ii)) & (pos_vec < bins(ii + 1));
    A_Maps(:,ii) = mean(proj(:,binMask), 2);  % since using mean no need for occ normalization. 
end

% normalize to peak in each map;
A_Maps = A_Maps./max(A_Maps,[],2);


[~, peak_idx] = max(A_Maps,[],2);
[~, s_idx] = sort(peak_idx); 

