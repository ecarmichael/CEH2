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
pos_Enc.label = {'x', 'y'}; 
pos_Enc.cfg.hdr{1} = behav_enc_a.json; 

pos_Rec = tsd(behav_rec_a.time, behav_rec_a.position');
pos_Rec.label = {'x', 'y'}; 
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
% Enc_ts = ts; 
% 
% for ii = size(bin_Enc_trl.data,1):-1:1
% Enc_ts.t{ii} = bin_Enc_trl.tvec(bin_Enc_trl.data(ii,:) ==1); 
% Enc_ts.label{ii} = num2str(ii); 
% end
% 
% cfg_pca = [];
% cfg_pca.plot = 1;
% cfg_pca.mov = 1;
% cfg_pca.bin_s = 0.5; 
% 
% [A_Temp, time_proj] = MS_PCA_ICA(cfg_pca, Enc_ts, pos_Enc_trl); 
% 
% 
% figure(999)
% for ii = 1:size(time_proj, 1)
%     
%   
%     plot(pos_Enc_trl.data(1,:), pos_Enc_trl.data(2,:), '.k')
%     hold
%         keep_idx = time_proj(1,:)> 8; 
%     plot(pos_Enc_trl.data(1,keep_idx), pos_Enc_trl.data(2,keep_idx), '.r')
%     title(num2str(ii))
%     
%     
%     
% end

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

lin_spd = getLinSpd([], pos_Enc_trl); 

cfg = []; cfg.method = 'raw'; cfg.operation = '>'; cfg.threshold = 2.5; % speed limit in cm/sec
iv_fast = TSDtoIV(cfg,lin_spd); % only keep intervals with speed above thresh

cfg = []; cfg.method = 'raw'; cfg.operation = '<'; cfg.threshold = 35; % speed limit in cm/sec
iv_max = TSDtoIV(cfg,lin_spd); % only keep intervals with speed above thresh

pos_enc_trl_mov = restrict(pos_Enc_trl, iv_fast);
pos_enc_trl_mov = restrict(pos_Enc_trl, iv_max);

bin_enc_trl_mov = restrict(bin_Enc_trl, iv_fast);
bin_enc_trl_mov = restrict(bin_Enc_trl, iv_max);


for kk = 1:19:size(bin_enc_trl_mov.data, 1)
c_idx = kk:kk+19;  

figure(kk)
for ii = 1:length(c_idx)
    subplot(4,5,ii)
    
    plot(pos_enc_trl_mov.data(1,:), pos_enc_trl_mov.data(2,:), '.k')
    hold
    keep_idx = nearest_idx3(bin_enc_trl_mov.tvec(logical(bin_enc_trl_mov.data(c_idx(ii),:))), pos_enc_trl_mov.tvec); 
    s = scatter(pos_enc_trl_mov.data(1,keep_idx), pos_enc_trl_mov.data(2,keep_idx), 20,pos_enc_trl_mov.tvec(keep_idx));
    s.MarkerFaceColor = s.MarkerEdgeColor; 
    title(num2str(c_idx(ii)))
    

end

end


%% rate maps
% set up bins
SET_xmin = floor(min(pos_enc_trl_mov.data(2,:))/100)*100 ; SET_ymin = floor(min(pos_enc_trl_mov.data(1,:))/100)*100; % set up bins
SET_xmax = ceil(max(pos_enc_trl_mov.data(2,:))/100)*100; SET_ymax = ceil(max(pos_enc_trl_mov.data(1,:))/100)*100;

% SET_xmin = 100; 
% SET_xmax = 400;
% SET_ymin = 200;
% SET_ymax = 700;
SET_xBinSz =5; SET_yBinSz = 5;


 
x_edges = SET_xmin:SET_xBinSz:SET_xmax;
y_edges = SET_ymin:SET_yBinSz:SET_ymax;

% gaussian kernal
kernel = gausskernel([1 1],1); % 2d gaussian in bins

% compute occupancy encode
occ_hist = histcn(pos_enc_trl_mov.data(1:2,:)',y_edges,x_edges); % 2-D version of histc()
occ_hist = conv2(occ_hist,kernel,'same');
 
no_occ_idx = find(occ_hist < 15); % NaN out bins never visited
occ_hist(no_occ_idx) = NaN;
 
occ_hist = occ_hist .* (mode(diff(pos_enc_trl_mov.tvec))); % convert samples to seconds using video frame rate (30 Hz)



m = 4;
n = 4;
s1_idx = 1:2:n*m; 
s2_idx = 2:2:n*m; 
figure(102)
clf
ip = 0; 
for ii = 1:size(bin_enc_trl_mov.data,1)
%     keep_idx = nearest_idx3(bin_enc_trl_mov.tvec(logical(bin_enc_trl_mov.data(c_idx(ii),:))), pos_enc_trl_mov.tvec); 

    spk_x = interp1(pos_enc_trl_mov.tvec,pos_enc_trl_mov.data(1,:),bin_enc_trl_mov.tvec(logical(bin_enc_trl_mov.data((ii),:))),'linear');
    spk_y = interp1(pos_enc_trl_mov.tvec,pos_enc_trl_mov.data(2,:),bin_enc_trl_mov.tvec(logical(bin_enc_trl_mov.data((ii),:))),'linear');
    if ip >= (n*m)/2
        figure(102+ii)
        ip = 1;
    else
        ip = ip+1;
    end
    
    disp(ip)
    subplot(m, n, s1_idx(ip))
    plot(pos_enc_trl_mov.data(1,:), pos_enc_trl_mov.data(2,:), '.k');
    title(num2str(ii))
    hold on
    plot(spk_x, spk_y, '.r');axis off;
    
    % compute the rate map
    if isempty(spk_x) || isempty(spk_y)
        continue
    end
    spk_hist = histcn([spk_x, spk_y],y_edges,x_edges);
    spk_hist = conv2(spk_hist,kernel, 'same');
    
    spk_hist(no_occ_idx) = NaN;
    tc = spk_hist./occ_hist;

    subplot(m, n, s2_idx(ip))
    pcolor(tc'); shading flat; axis off; %colorbar('Location', 'northoutside')
   
end

%% save for later. 