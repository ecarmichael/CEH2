function [sig_act] = MS_Asmbly_activations(proj, threshold, min_dist)




%% check for sig reactivations
sig_act = [];
for ii = size(proj_out,1):-1:1
    [~, p_idx] = findpeaks((proj(ii,:)),'MinPeakHeight', threshold ,'MinPeakDistance', min_dist);
    
    if ~isempty(p_idx)
        sig_act(ii) = length(p_idx);
    else
        sig_act(ii) = NaN;
    end
    
end