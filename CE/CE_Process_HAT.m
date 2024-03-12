function CE_Process_HAT(Ca_dir, nlx_dir)
%% CE_Process_HAT: loads, aligns, and runs sleep analyses on v4 HAT data



%% load the NLX data

cd(nlx_dir)
Meta = MS_Load_meta; 

cfg_csc = [];
cfg_csc.fc = {Meta.goodCSC, Meta.goodCSC2}; 
cfg_csc.desired_sampling_frequency = 2000; 

csc = MS_LoadCSC(cfg_csc);


evts = MS_LoadEvents();

% find the Sync TTLs

r_start_idx = find(contains(evts.label, 'Starting Recording')); 

r_stop_idx = find(contains(evts.label, 'Stopping Recording')); 

sync_ttl = sort(unique([evts.t{3}, evts.t{4}]));

start_idx = [1 find(diff(sync_ttl) > 5)+1]; 
end_idx = [find(diff(sync_ttl) > 5) length(sync_ttl)]; 


csc_pre = restrict(csc, sync_ttl(start_idx(1)), sync_ttl(end_idx(1))); 

csc_trk = restrict(csc, sync_ttl(start_idx(2)), sync_ttl(end_idx(2))); 

csc_post = restrict(csc, sync_ttl(start_idx(3)), sync_ttl(end_idx(3))); 

ttl_pre = restrict(evts, sync_ttl(1), 

%% load the concatenated ms file


load([Ca_dir filesep 'ms.mat'])

length(ms.RawTraces)




