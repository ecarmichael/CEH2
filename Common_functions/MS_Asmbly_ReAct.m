function [proj_out, data_h, tvec, R_threshold] = MS_Asmbly_ReAct(data_in, Temp_in,ms, bin_size);



this_time =  0:1/mode(diff(ms.time)):(length(data_in)/mode(diff(ms.time)));
this_time = this_time(1:end-1);

tbin_edges = 0:bin_size:(length(data_in)/mode(diff(ms.time))); % vector of time bin edges (for histogram)
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

%%
rng(123,'twister')
nShuff = 500;
shuff_mat = [];

for iS = 1:nShuff
    tic
    shuff_data = NaN(size(data_h));
    for ic = 1:size(data_h,2)
        
        vals = (size(data_h,1)-1).*rand(1,1) + 1;

        shuff_data(:,ic) = circshift(data_h(:,ic), floor(vals));
    end
    
    shuff_proj = assembly_activity(Temp_in,shuff_data');
    
    shuff_mat = [shuff_mat; shuff_proj];
    
   
    
end

R_threshold = prctile(shuff_mat(shuff_mat >0), 99, 'all');

