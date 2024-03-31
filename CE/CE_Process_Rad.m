function CE_Process_Rad(Ca_dir, nlx_dir)
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


csc_enc = restrict(csc, sync_ttl(start_idx(1)), sync_ttl(end_idx(1))); 
evt_enc = restrict(evts, sync_ttl(start_idx(1)), sync_ttl(end_idx(1))); 
ttl_enc = sort(unique([evt_enc.t{3}, evt_enc.t{4}]));

csc_sleep = restrict(csc, sync_ttl(start_idx(2)), sync_ttl(end_idx(2))); 
evt_sleep = restrict(evts, sync_ttl(start_idx(2)), sync_ttl(end_idx(2))); 
ttl_sleep = sort(unique([evt_sleep.t{3}, evt_sleep.t{4}]));

csc_rec = restrict(csc, sync_ttl(start_idx(3)), sync_ttl(end_idx(3))); 
evt_rec = restrict(evts, sync_ttl(start_idx(3)), sync_ttl(end_idx(3))); 
ttl_rec = sort(unique([evt_rec.t{3}, evt_rec.t{4}]));


%% load the concatenated ms file


load([Ca_dir filesep 'ms.mat'])


%% align the nlx tvecs

if length(ttl_enc) - length(ms.tvecs{1}) > 0
    fprintf('<strong>%s: NLX: %.0f   -    %.0f (miniscope) | diff: %.0f</strong>', 'pre', length(ttl_enc), length(ms.tvecs{1}), length(ttl_enc) - length(ms.tvecs{1}))
    
    nlx_temp = ttl_enc - ttl_enc(1);
    
    ttl_enc = interp1(nlx_temp, ttl_enc, ms.tvecs{1} - ms.tvecs{1}(1));
    
    fprintf('<strong>... fixed --> diff = %.0f</strong>\n', length(ttl_enc) - length(ms.tvecs{1}))
else
    fprintf('%s: NLX: %.0f   -    %.0f (miniscope) | diff: %.0f\n', 'pre', length(ttl_enc), length(ms.tvecs{1}), length(ttl_enc) - length(ms.tvecs{1}))
end

if length(ttl_sleep) - length(ms.tvecs{2}) > 0
    fprintf('<strong>%s: NLX: %.0f   -    %.0f (miniscope) | diff: %.0f</strong>', 'trk', length(ttl_sleep), length(ms.tvecs{2}), length(ttl_sleep) - length(ms.tvecs{2}))
    
    nlx_temp = ttl_sleep - ttl_sleep(1);
    
    ttl_sleep = interp1(nlx_temp, ttl_sleep, ms.tvecs{2} - ms.tvecs{2}(1));
    
    fprintf('<strong>... fixed --> diff = %.0f</strong>\n', length(ttl_sleep) - length(ms.tvecs{2}))
else
    fprintf('%s: NLX: %.0f   -    %.0f (miniscope) | diff: %.0f\n', 'trk', length(ttl_sleep), length(ms.tvecs{2}), length(ttl_sleep) - length(ms.tvecs{2}))
end


if length(ttl_rec) - length(ms.tvecs{3}) > 0
    fprintf('%s: NLX: %.0f   -    %.0f (miniscope) | diff: %.0f', 'post', length(ttl_rec), length(ms.tvecs{3}), length(ttl_rec) - length(ms.tvecs{3}))
    
    nlx_temp = ttl_rec - ttl_rec(1);
    
    ttl_rec = interp1(nlx_temp, ttl_rec, ms.tvecs{3} - ms.tvecs{3}(1));
    
    fprintf('<strong>... fixed --> diff = %.0f</strong>\n', length(ttl_rec) - length(ms.tvecs{3}))
else
    fprintf('%s: NLX: %.0f   -    %.0f (miniscope) | diff: %.0f\n', 'post', length(ttl_rec), length(ms.tvecs{3}), length(ttl_rec) - length(ms.tvecs{3}))
end


%% align the miniscope motion to the LFP


sleep_motion = interp1(ms.HD{2}.tvec - ms.HD{2}.tvec(1), ms.HD{2}.motion, csc_sleep.tvec - csc_sleep.tvec(1))'; 


%% get the hypnograms


hypno = MS_get_hypno(csc_sleep, sleep_motion,[], 75); 
saveas(gcf, [Ca_dir filesep 'hypno.fig']); 
saveas(gcf, [Ca_dir filesep 'hypno.png']); 
save('hypno.mat', 'hypno', '-v7.3')
close(gcf); 



clear sleep_motion 

%% cut the ms struct into REM SWS
temp_raw = ms.RawTraces(1:ms.timestamps(1)+1,:);
temp_bin = ms.Binary(1:ms.timestamps(1)+1,:);

REM.RawTraces = []; REM.tvec = []; 
REM.Binary = []; 
for ii =  1:length(hypno.cfg.REM.tstart)
    s_idx = nearest_idx(hypno.cfg.REM.tstart(ii)-csc_enc.tvec(1), ms.tvecs{1}); 
    e_idx = nearest_idx(hypno.cfg.REM.tend(ii)-csc_enc.tvec(1), ms.tvecs{1}); 

    REM.RawTraces = [REM.RawTraces;  temp_raw(s_idx:e_idx,:)];
    REM.Binary = [REM.Binary;  temp_bin(s_idx:e_idx,:)];

    REM.tvec = [REM.tvec; ms.tvecs{1}(s_idx:e_idx)];
    
    REM.split_idx(ii) = length(s_idx:e_idx); 

end



SWS.RawTraces = []; SWS.tvec = []; 
SWS.Binary = []; 
for ii =  1:length(hypno.cfg.SWS.tstart)
    s_idx = nearest_idx(hypno.cfg.SWS.tstart(ii)-csc_enc.tvec(1), ms.tvecs{1}); 
    e_idx = nearest_idx(hypno.cfg.SWS.tend(ii)-csc_enc.tvec(1), ms.tvecs{1}); 

    SWS.RawTraces = [SWS.RawTraces;  temp_raw(s_idx:e_idx,:)];
    SWS.Binary = [SWS.Binary;  temp_bin(s_idx:e_idx,:)];

    SWS.tvec = [SWS.tvec; ms.tvecs{1}(s_idx:e_idx)];
    
    SWS.split_idx(ii) = length(s_idx:e_idx); 

end



%% split the awake phase into trials
rad_name = dir('Radial_log*.m'); 

run(rad_name.name)

target = Rad.(['D' rad_name.name(17:end-2)]).correct;
trl = Rad.(['D' rad_name.name(17:end-2)]).(['m' rad_name.name(12:15)]);

% Encoding trials

Enc_iv = iv(trl.encode.tstart+sync_ttl(start_idx(1)), trl.encode.tend+sync_ttl(start_idx(1)));
Enc_iti_iv = iv(trl.encode.tstart-60+sync_ttl(start_idx(1)), trl.encode.tstart+sync_ttl(start_idx(1)));


csc_enc_trl = restrict(csc_enc, Enc_iv);
enc_tstart = csc_enc_trl.tvec(1); 
csc_enc_trl.tvec = csc_enc_trl.tvec - enc_tstart; 

csc_enc_iti = restrict(csc_enc, Enc_iti_iv);
csc_enc_iti.tvec = csc_enc_iti.tvec - enc_tstart; 






