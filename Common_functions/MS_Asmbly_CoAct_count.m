function [alt_sig, base_sig] = MS_Asmbly_CoAct_count(base_mat, alt_mat, thresh)
%% 

if nargin < 3
    thresh = 1.96;
end

sig_idx =  base_mat > thresh; 
alt_sig = alt_mat > thresh;


    alt_sig = (sum(alt_sig, 'all') - numel(diag(sig_idx))); 
    
    
    base_sig = (sum(sig_idx, 'all') - numel(diag(sig_idx)));
    
%     fprintf('Alt_mat had %.0f (%.0f%%) significant assemblies of the %0.0f found in base_mat \n',...
%     alt_sig, ((sum(alt_sig, 'all') - numel(diag(sig_idx))) / (numel(sig_idx) - numel(diag(sig_idx))))*100,...
%     base_sig)


