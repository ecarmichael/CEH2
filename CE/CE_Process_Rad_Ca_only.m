function CE_Process_Rad_Ca_only(Ca_dir)
%% CE_Process_HAT: loads, aligns, and runs sleep analyses on v4 HAT data

% Note:  conv_fact for DLC should be [9.43 9.43]; 

%% load the concatenated ms file
warning off
tic
load([Ca_dir filesep 'ms.mat'], 'ms')
toc
warning on

if exist([Ca_dir    filesep 'keep_idx.mat'], 'file')
    load([Ca_dir    filesep 'keep_idx.mat'], 'keep_idx')
    cfg_ms = [];
    cfg_ms.remove_idx = ~keep_idx;
    ms = MS_Remove_trace(cfg_ms, ms);

end


% split into encoding and recall

ms_enc = MS_restrict(ms, ms.time(1), ms.time(ms.timestamps(1)+1)); 
ms_rec = MS_restrict(ms, ms.time(ms.timestamps(1)+2), ms.time(end)); 

ms_rec_tstart = ms_rec.time(1); 
ms_rec.time = ms_rec.time - ms_rec.time(1); 

%% load the behaviour

load([Ca_dir filesep 'behav_enc.mat'])
load([Ca_dir filesep 'behav_rec.mat'])


%% align the behaviour to the Ca time

% encoding
behav_enc_a = MS_align_data(behav_enc, ms_enc); 
behav_rec_a = MS_align_data(behav_rec, ms_rec); 



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