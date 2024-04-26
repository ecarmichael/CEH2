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

ms_Enc = MS_restrict(ms, ms.time(1), ms.time(ms.timestamps(1)+1)); 
ms_Rec = MS_restrict(ms, ms.time(ms.timestamps(1)+2), ms.time(end)); 

ms_rec_tstart = ms_Rec.time(1); 
ms_Rec.time = ms_Rec.time - ms_Rec.time(1); 

%% load the behaviour

load([Ca_dir filesep 'behav_enc.mat'])
load([Ca_dir filesep 'behav_rec.mat'])


%% align the behaviour to the Ca time

% encoding
behav_enc_a = MS_align_data(behav_enc, ms_Enc); 
behav_rec_a = MS_align_data(behav_rec, ms_Rec); 


%% convert to tsd for restricting data

pos_Enc = tsd(behav_enc_a.time, behav_enc_a.position');  
pos_Enc.cfg.hdr{1} = behav_enc_a.json; 

pos_Rec = tsd(behav_rec_a.time, behav_rec_a.position');
pos_Rec.cfg.hdr{1} = behav_rec_a.json; 

ms_Enc = tsd(ms_Enc.time, ms_Enc.Binary');
ms_Enc.cfg.hdr{1} = ms.Exp_json{1}; 

ms_Rec = tsd(ms_Rec.time, ms_Rec.Binary');
ms_Rec.cfg.hdr{1} = ms.Exp_json{2}; 



%% split the awake phase into trials
rad_name = dir([Ca_dir filesep 'Radial_log*.m']); 

run([rad_name.folder filesep  rad_name.name])

target = Rad.(['D' rad_name.name(17:end-2)]).correct;
trl = Rad.(['D' rad_name.name(17:end-2)]).(['m' rad_name.name(12:15)]);

% Encoding trials

Enc_iv = iv(trl.encode.tstart(1:4), trl.encode.tend(1:4));

iti_s = [trl.encode.tstart(1:4)-60, trl.encode.tstart(5)]; 
iti_e = [trl.encode.tstart(1:4), trl.encode.tend(5)]; 

Enc_iti_iv = iv(iti_s, iti_e);

% Recall trials

Rec_iv = iv(trl.recall.tstart(1:4), trl.recall.tend(1:4));

iti_s = [trl.recall.tstart(1:4)-60, trl.recall.tstart(5)]; 
iti_e = [trl.recall.tstart(1:4), trl.recall.tend(5)]; 

Rec_iti_iv = iv(iti_s, iti_e);


%% restrict the behaviour to trials

pos_Enc_trl = restrict(pos_Enc, Enc_iv);

pos_Enc_iti = restrict(pos_Enc, Enc_iti_iv);

pos_Rec_trl = restrict(pos_Rec, Rec_iv);

pos_Rec_iti = restrict(pos_Rec, Rec_iti_iv);

%% restrict the calcium data to the trial time

% bin_enc = tsd(ms.tvecs{1}+sync_ttl(start_idx(1)), ms.Binary(1:length(ms.tvecs{1}),:)'); 
% bin_enc.cfg.hdr{1} = ms.Exp_json{1}; 
% 
% bin_rec = tsd(ms.tvecs{3}+sync_ttl(start_idx(3)), ms.Binary((length(ms.tvecs{1})+length(ms.tvecs{2}))+1:end,:)'); 
% bin_rec.cfg.hdr{1} = ms.Exp_json{3}; 



bin_Enc_trl = restrict(ms_Enc, Enc_iv);

bin_Enc_iti = restrict(ms_Enc, Enc_iti_iv);


bin_Rec_trl = restrict(ms_Rec, Rec_iv);

bin_Rec_iti = restrict(ms_Rec, Rec_iti_iv);



%% assembly dectection. screening
Enc_ts = ts; 

for ii = size(bin_Enc_trl.data,1):-1:1
Enc_ts.t{ii} = bin_Enc_trl.tvec(bin_Enc_trl.data(ii,:) ==1); 
Enc_ts.label{ii} = num2str(ii); 
end

cfg_pca = [];
cfg_pca.plot = 1;
cfg_pca.mov = 1;
cfg_pca.bin_s = 0.5; 

[A_Temp, time_proj] = MS_PCA_ICA(cfg_pca, Enc_ts, pos_Enc_trl); 


figure(999)
for ii = 1:size(time_proj, 1)
    
  
    plot(pos_Enc_trl.data(1,:), pos_Enc_trl.data(2,:), '.k')
    hold
        keep_idx = time_proj(1,:)> 8; 
    plot(pos_Enc_trl.data(1,keep_idx), pos_Enc_trl.data(2,keep_idx), '.r')
    title(num2str(ii))
    
    
    
end

%% check aligment

figure(191)
clf
ax(1) = subplot(4,2,1); 
plot(pos_Enc_trl.tvec, pos_Enc_trl.data(1,:), 'b')
hold on
plot(pos_Enc_iti.tvec, pos_Enc_iti.data(1,:), 'r')


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

c_idx = 110:120; 

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


%% save for later. 