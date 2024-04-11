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
warning off
tic
load([Ca_dir filesep 'ms.mat'], 'ms')
toc
warning on

%% load the behaviour

load([Ca_dir filesep 'behav_enc.mat'])
load([Ca_dir filesep 'behav_rec.mat'])

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
temp_raw = ms.RawTraces(1:ms.timestamps(2)+1,:);
temp_bin = ms.Binary(1:ms.timestamps(2)+1,:);

REM.RawTraces = []; REM.tvec = []; 
REM.Binary = []; 
for ii =  1:length(hypno.cfg.REM.tstart)
    s_idx = nearest_idx3(hypno.cfg.REM.tstart(ii)-csc_sleep.tvec(1), ms.tvecs{2} - ms.tvecs{2}(1)); 
    e_idx = nearest_idx3(hypno.cfg.REM.tend(ii)-csc_sleep.tvec(1), ms.tvecs{2} - ms.tvecs{2}(1)); 

    REM.RawTraces = [REM.RawTraces;  temp_raw(s_idx:e_idx,:)];
    REM.Binary = [REM.Binary;  temp_bin(s_idx:e_idx,:)];

    REM.tvec = [REM.tvec; ms.tvecs{2}(s_idx:e_idx)];
    
    REM.split_idx(ii) = length(s_idx:e_idx); 

end



SWS.RawTraces = []; SWS.tvec = []; 
SWS.Binary = []; 
for ii =  1:length(hypno.cfg.SWS.tstart)
    s_idx = nearest_idx(hypno.cfg.SWS.tstart(ii)-csc_sleep.tvec(1), ms.tvecs{2}- ms.tvecs{2}(1)); 
    e_idx = nearest_idx(hypno.cfg.SWS.tend(ii)-csc_sleep.tvec(1), ms.tvecs{2}- ms.tvecs{2}(1)); 

    SWS.RawTraces = [SWS.RawTraces;  temp_raw(s_idx:e_idx,:)];
    SWS.Binary = [SWS.Binary;  temp_bin(s_idx:e_idx,:)];

    SWS.tvec = [SWS.tvec; ms.tvecs{2}(s_idx:e_idx)];
    
    SWS.split_idx(ii) = length(s_idx:e_idx); 

end



%% split the awake phase into trials
rad_name = dir([Ca_dir filesep 'Radial_log*.m']); 

run([rad_name.folder filesep  rad_name.name])

target = Rad.(['D' rad_name.name(17:end-2)]).correct;
trl = Rad.(['D' rad_name.name(17:end-2)]).(['m' rad_name.name(12:15)]);

% Encoding trials

Enc_iv = iv(trl.encode.tstart(1:4)+sync_ttl(start_idx(1)), trl.encode.tend(1:4)+sync_ttl(start_idx(1)));

iti_s = [trl.encode.tstart(1:4)-60+sync_ttl(start_idx(1)), trl.encode.tstart(5)+sync_ttl(start_idx(1))]; 
iti_e = [trl.encode.tstart(1:4)+sync_ttl(start_idx(1)), trl.encode.tend(5)+sync_ttl(start_idx(1))]; 

Enc_iti_iv = iv(iti_s, iti_e);

% Recall trials

Rec_iv = iv(trl.recall.tstart(1:4)+sync_ttl(start_idx(3)), trl.recall.tend(1:4)+sync_ttl(start_idx(3)));

iti_s = [trl.recall.tstart(1:4)-60+sync_ttl(start_idx(3)), trl.recall.tstart(5)+sync_ttl(start_idx(3))]; 
iti_e = [trl.recall.tstart(1:4)+sync_ttl(start_idx(3)), trl.recall.tend(5)+sync_ttl(start_idx(3))]; 

Rec_iti_iv = iv(iti_s, iti_e);

enc_t0 = sync_ttl(start_idx(1));
rec_t0 = sync_ttl(start_idx(3));

%% restrict the lfp data

% csc encoding
csc_Enc_trl = restrict(csc_enc, Enc_iv);
csc_Enc_trl.tvec = csc_Enc_trl.tvec - enc_t0; 

%enc iti
csc_Enc_iti = restrict(csc_enc, Enc_iti_iv);
csc_Enc_iti.tvec = csc_Enc_iti.tvec - enc_t0; 


% csc recall
csc_Rec_trl = restrict(csc_rec, Rec_iv);
csc_Rec_trl.tvec = csc_Rec_trl.tvec - rec_t0; 

%rec iti
csc_Rec_iti = restrict(csc_rec, Rec_iti_iv);
csc_Rec_iti.tvec = csc_Rec_iti.tvec - rec_t0; 


% figure(1010)
% subplot(2,1,1)
% pl_cfg = [];
% PlotTSDfromIV(pl_cfg, Enc_iv, csc_enc)
% ylabel('Trials')
% 
% subplot(2,1,2)
% pl_cfg = [];
% PlotTSDfromIV(pl_cfg, Enc_iti_iv, csc_enc)
% ylabel('Trials')
%% restrict the behaviour to trials

pos_enc = tsd(behav_enc.time+sync_ttl(start_idx(1)), [behav_enc.position(:,1), behav_enc.position(:,2)]'); 
pos_enc.cfg.hdr{1} = behav_enc.json; 

pos_rec = tsd(behav_rec.time+sync_ttl(start_idx(3)), [behav_rec.position(:,1), behav_rec.position(:,2)]'); 
pos_rec.cfg.hdr{1} = behav_rec.json; 


pos_Enc_trl = restrict(pos_enc, Enc_iv);
pos_Enc_trl.tvec = pos_Enc_trl.tvec - enc_t0; 

pos_Enc_iti = restrict(pos_enc, Enc_iti_iv);
pos_Enc_iti.tvec = pos_Enc_iti.tvec - enc_t0; 

pos_Rec_trl = restrict(pos_rec, Rec_iv);
pos_Rec_trl.tvec = pos_Rec_trl.tvec - rec_t0; 

pos_Rec_iti = restrict(pos_rec, Rec_iti_iv);
pos_Rec_iti.tvec = pos_Rec_iti.tvec - rec_t0; 

%% restrict the calcium data to the trial time

% bin_enc = tsd(ms.tvecs{1}+sync_ttl(start_idx(1)), ms.Binary(1:length(ms.tvecs{1}),:)'); 
% bin_enc.cfg.hdr{1} = ms.Exp_json{1}; 
% 
% bin_rec = tsd(ms.tvecs{3}+sync_ttl(start_idx(3)), ms.Binary((length(ms.tvecs{1})+length(ms.tvecs{2}))+1:end,:)'); 
% bin_rec.cfg.hdr{1} = ms.Exp_json{3}; 

bin_tsd =  tsd(ms.time+sync_ttl(start_idx(1)), ms.Binary'); 


bin_Enc_trl = restrict(bin_tsd, Enc_iv);
bin_Enc_trl.tvec = bin_Enc_trl.tvec - enc_t0; 

bin_Enc_iti = restrict(bin_tsd, Enc_iti_iv);
bin_Enc_iti.tvec = bin_Enc_iti.tvec - enc_t0; 


bin_Rec_trl = restrict(bin_tsd, Rec_iv);
bin_Rec_trl.tvec = bin_Rec_trl.tvec - rec_t0; 

bin_Rec_iti = restrict(bin_tsd, Rec_iti_iv);
bin_Rec_iti.tvec = bin_Rec_iti.tvec - rec_t0; 


%% check aligment

figure(191)
clf
ax(1) = subplot(4,2,1); 
plot(pos_Enc_trl.tvec, pos_Enc_trl.data(1,:), 'b')
hold on
plot(pos_Enc_iti.tvec, pos_Enc_iti.data(1,:), 'r')

ax(2) = subplot(4,2,3); 
plot(csc_Enc_trl.tvec, csc_Enc_trl.data(1,:), 'b')
hold on
plot(csc_Enc_iti.tvec, csc_Enc_iti.data(1,:), 'r')

ax(3) = subplot(4,2,5);
scatter(pos_Enc_trl.tvec, pos_Enc_trl.data(1,:), 'b')
hold on
scatter(pos_Enc_iti.tvec, pos_Enc_iti.data(1,:), 'r')

ax(4) = subplot(4,2,7);
plot(bin_Enc_trl.tvec, bin_Enc_trl.data(100,:), 'b')
hold on
plot(bin_Enc_iti.tvec, bin_Enc_iti.data(100,:), 'r')

linkaxes(ax, 'x'); 
xlim([bin_Enc_iti.tvec(1) bin_Enc_iti.tvec(end)])

% recall
ax2(1) = subplot(4,2,2); 
plot(pos_Rec_trl.tvec, pos_Rec_trl.data(1,:), 'b')
hold on
plot(pos_Rec_iti.tvec, pos_Rec_iti.data(1,:), 'r')

ax2(2) = subplot(4,2,4); 
plot(csc_Rec_trl.tvec, csc_Rec_trl.data(1,:), 'b')
hold on
plot(csc_Rec_iti.tvec, csc_Rec_iti.data(1,:), 'r')

ax2(3) = subplot(4,2,6);
scatter(pos_Rec_trl.tvec, pos_Rec_trl.data(1,:), 'b')
hold on
scatter(pos_Rec_iti.tvec, pos_Rec_iti.data(1,:), 'r')

ax2(4) = subplot(4,2,8);
plot(bin_Rec_trl.tvec, bin_Rec_trl.data(100,:), 'b')
hold on
plot(bin_Rec_iti.tvec, bin_Rec_iti.data(100,:), 'r')

linkaxes(ax2, 'x'); 
xlim([bin_Rec_iti.tvec(1) bin_Rec_iti.tvec(end)])



%% check for place cells
figure(199)
clf

c_idx = 230:246; 

for ii = 1:length(c_idx)
    subplot(4,4,ii)
    
    plot(pos_Enc_trl.data(1,:), pos_Enc_trl.data(2,:), '.k')
    hold
    keep_idx = nearest_idx3(bin_Enc_trl.tvec(logical(bin_Enc_trl.data(c_idx(ii),:))), pos_Enc_trl.tvec); 
    plot(pos_Enc_trl.data(1,keep_idx), pos_Enc_trl.data(2,keep_idx), '.r')
    title(num2str(c_idx(ii)))
    

end



for ii = 1:length(c_idx)
    
    
    
    
end

%% try to plot the power of different lfp elemets

figure(99)
subplot

%% save for later. 