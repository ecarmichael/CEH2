function [proj_out,ReAct_stats, data_h, tvec, shuff] = MS_Asmbly_ReAct(cfg, data_in, Temp_in,ms, bin_size);




%% initialize
if isempty(cfg)
    cfg.nShuff = 500;
    cfg.thresh = 99; % in percentile
    cfg.ms_fs = []; 
end

if ~isfield(cfg, 'ms_fs') || isempty(cfg.ms_fs)
    
    ms_fs = mode(diff(ms.time)); 
    
else
    
ms_fs = cfg.ms_fs; 
    
end

%% coner the data into binned data

this_time =  0:1/ms_fs:(length(data_in)/ms_fs);
this_time = this_time(1:end-1);

tbin_edges = 0:bin_size:(length(data_in)/ms_fs); % vector of time bin edges (for histogram)
tbin_centers = tbin_edges(1:end-1)+bin_size/2; % vector of time bin centers (for plotting)


data_h= [];
for ii = size(data_in,2):-1:1
    

        this_cell = this_time(find(data_in(:,ii)));
        
    spk_count = histc(this_cell,tbin_edges); % get spike counts for each bin
    spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.
%     data_h_rem(:,ii) = zscore(spk_count);
    data_h(:,ii) = spk_count;

end
tvec = tbin_centers;

%% get the reactivation projection
rng(123, 'twister')

proj_out = assembly_activity(Temp_in ,data_h');


%% use a shuffle permutation to get the R_threshold and null 

[ReAct_stats, shuff.data, shuff.proj] = MS_Asmbly_proj_thresh(data_h, Temp_in, cfg.nShuff, cfg.thresh); 


%% get the reactivation stats
    
ReAct_stats.p_val = [];
ReAct_stats.rate = [];
ReAct_stats.rate_p = [];



for ii = size(proj_out,1):-1:1
    ReAct_stats.p_val(ii) = sum(sum(shuff.data > ReAct_stats.R_thresh,2) > sum(proj_out(ii,:) > ReAct_stats.R_thresh))/ size(shuff.data,1);
    
    ReAct_stats.rate(ii) = sum(proj_out(ii,:) > ReAct_stats.R_thresh) / ((tvec(end) - tvec(1))/60);
    
    ReAct_stats.shuff_rate = sum(shuff.data > ReAct_stats.R_thresh,2)./ ((tvec(end) - tvec(1))/60);
    
    ReAct_stats.rate_p(ii) = sum(ReAct_stats.shuff_rate > ReAct_stats.rate(ii)) / length(ReAct_stats.shuff_rate);
    
end
    
