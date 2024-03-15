function CE_Process_HAT(Ca_dir, nlx_dir)
%% CE_Process_HAT: loads, aligns, and runs sleep analyses on v4 HAT data


%% load the NLX data

cd(nlx_dir)
Meta = MS_Load_meta; 

cfg_csc = [];
cfg_csc.fc = {Meta.goodCSC}; 
cfg_csc.desired_sampling_frequency = 2000; 

csc = MS_LoadCSC(cfg_csc);


evts = MS_LoadEvents();
evts.cfg.history.mfun = [];
evts.cfg.history.cfg = []; 
% find the Sync TTLs

r_start_idx = find(contains(evts.label, 'Starting Recording')); 

r_stop_idx = find(contains(evts.label, 'Stopping Recording')); 

sync_ttl = sort(unique([evts.t{3}, evts.t{4}]));

start_idx = [1 find(diff(sync_ttl) > 5)+1]; 
end_idx = [find(diff(sync_ttl) > 5) length(sync_ttl)]; 


csc_pre = restrict(csc, sync_ttl(start_idx(1)), sync_ttl(end_idx(1))); 
evt_pre = restrict(evts, sync_ttl(start_idx(1)), sync_ttl(end_idx(1))); 
ttl_pre = sort(unique([evt_pre.t{3}, evt_pre.t{4}]));

csc_trk = restrict(csc, sync_ttl(start_idx(2)), sync_ttl(end_idx(2))); 
evt_trk = restrict(evts, sync_ttl(start_idx(2)), sync_ttl(end_idx(2))); 
ttl_trk = sort(unique([evt_trk.t{3}, evt_trk.t{4}]));

csc_post = restrict(csc, sync_ttl(start_idx(3)), sync_ttl(end_idx(3))); 
evt_post = restrict(evts, sync_ttl(start_idx(3)), sync_ttl(end_idx(3))); 
ttl_post = sort(unique([evt_post.t{3}, evt_post.t{4}]));


%% load the concatenated ms file


load([Ca_dir filesep 'ms.mat'])


%% align the nlx tvecs

if length(ttl_pre) - length(ms.tvecs{1}) > 0
    fprintf('<strong>%s: NLX: %.0f   -    %.0f (miniscope) | diff: %.0f</strong>', 'pre', length(ttl_pre), length(ms.tvecs{1}), length(ttl_pre) - length(ms.tvecs{1}))
    
    nlx_temp = ttl_pre - ttl_pre(1);
    
    ttl_pre = interp1(nlx_temp, ttl_pre, ms.tvecs{1} - ms.tvecs{1}(1));
    
    fprintf('<strong>... fixed --> diff = %.0f</strong>\n', length(ttl_pre) - length(ms.tvecs{1}))
else
    fprintf('%s: NLX: %.0f   -    %.0f (miniscope) | diff: %.0f\n', 'pre', length(ttl_pre), length(ms.tvecs{1}), length(ttl_pre) - length(ms.tvecs{1}))
end

if length(ttl_trk) - length(ms.tvecs{2}) > 0
    fprintf('<strong>%s: NLX: %.0f   -    %.0f (miniscope) | diff: %.0f</strong>', 'trk', length(ttl_trk), length(ms.tvecs{2}), length(ttl_trk) - length(ms.tvecs{2}))
    
    nlx_temp = ttl_trk - ttl_trk(1);
    
    ttl_trk = interp1(nlx_temp, ttl_trk, ms.tvecs{2} - ms.tvecs{2}(1));
    
    fprintf('<strong>... fixed --> diff = %.0f</strong>\n', length(ttl_trk) - length(ms.tvecs{2}))
else
    fprintf('%s: NLX: %.0f   -    %.0f (miniscope) | diff: %.0f\n', 'trk', length(ttl_trk), length(ms.tvecs{2}), length(ttl_trk) - length(ms.tvecs{2}))
end


if length(ttl_post) - length(ms.tvecs{3}) > 0
    fprintf('%s: NLX: %.0f   -    %.0f (miniscope) | diff: %.0f\n', 'post', length(ttl_post), length(ms.tvecs{3}), length(ttl_post) - length(ms.tvecs{3}))
    
    nlx_temp = ttl_post - ttl_post(1);
    
    ttl_post = interp1(nlx_temp, ttl_post, ms.tvecs{3} - ms.tvecs{3}(1));
    
    fprintf('<strong>... fixed --> diff = %.0f</strong>\n', length(ttl_post) - length(ms.tvecs{3}))
else
    fprintf('%s: NLX: %.0f   -    %.0f (miniscope) | diff: %.0f\n', 'post', length(ttl_post), length(ms.tvecs{3}), length(ttl_post) - length(ms.tvecs{3}))
end


%% align the miniscope motion to the LFP


pre_motion = interp1(ms.HD{1}.tvec - ms.HD{1}.tvec(1), ms.HD{1}.motion, csc_pre.tvec - csc_pre.tvec(1))'; 

post_motion = interp1(ms.HD{3}.tvec - ms.HD{3}.tvec(1), ms.HD{3}.motion, csc_post.tvec - csc_post.tvec(1))'; 


%% get the hypnograms


pre_hypno = MS_get_hypno(csc_pre, pre_motion); 
saveas(gcf, [Ca_dir filesep 'pre_hypno.fig']); 
saveas(gcf, [Ca_dir filesep 'pre_hypno.png']); 
close(gcf); 

post_hypno = MS_get_hypno(csc_post, post_motion); 
saveas(gcf, [Ca_dir filesep 'post_hypno.fig']); 
saveas(gcf, [Ca_dir filesep 'post_hypno.png']); 
close(gcf)

clear pre_motion post_motion

%% cut the ms struct into REM SWS
temp_raw = ms.RawTraces(1:ms.timestamps(1)+1,:);
temp_bin = ms.Binary(1:ms.timestamps(1)+1,:);

REM_pre.RawTraces = []; REM_pre.tvec = []; 
REM_pre.Binary = []; 
for ii =  1:length(pre_hypno.cfg.REM.tstart)
    s_idx = nearest_idx(pre_hypno.cfg.REM.tstart(ii)-csc_pre.tvec(1), ms.tvecs{1}); 
    e_idx = nearest_idx(pre_hypno.cfg.REM.tend(ii)-csc_pre.tvec(1), ms.tvecs{1}); 

    REM_pre.RawTraces = [REM_pre.RawTraces;  temp_raw(s_idx:e_idx,:)];
    REM_pre.Binary = [REM_pre.Binary;  temp_bin(s_idx:e_idx,:)];

    REM_pre.tvec = [REM_pre.tvec; ms.tvecs{1}(s_idx:e_idx)];
    
    split_idx(ii) = length(s_idx:e_idx); 

end


