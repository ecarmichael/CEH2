%% sandbox_PCA_ICA_examples


addpath('/home/williamslab/Documents/Github/Dos-Santos Assembly ICA/Dos-Santos Assembly ICA')
%% %%%%% EPHYS %%%%%%%%%%%%
load('/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Radial/Inter/M30_2023_06_10_Rad5.mat')

%% limit to movement in awake
cfg = []; cfg.method = 'raw'; cfg.operation = '>'; cfg.threshold = 2.5; % speed limit in cm/sec
iv_fast = TSDtoIV(cfg,out.Encode.speed); % only keep intervals with speed above thresh

out.Encode.S = restrict(out.Encode.S, iv_fast);
out.Encode.pos = restrict(out.Encode.pos, iv_fast);

%% get some assembly patterns


cfg_pca = [];
cfg_pca.plot = 1; 
cfg_pca.pos_ass = 1;
cfg_pca.bin_s = .005;
[A_temp, T_proj] = MS_PCA_ICA(cfg_pca, out.Encode.S, out.Encode.pos, out.Encode.csc); 


%% use Tingley Peyrache code



cfg.bin_s = .01; 
    tstart = out.Encode.csc.tvec(1);
    tend = out.Encode.csc.tvec(end);  

tbin_edges_enc  = tstart:cfg.bin_s:tend; % vector of time bin edges (for histogram)
tbin_centers_enc = tbin_edges_enc (1:end-1)+cfg.bin_s/2; % vector of time bin centers (for plotting)

data_h = [];
for ii = size(out.Encode.S.t,2):-1:1
    
    
    spk_count = histc(out.Encode.S.t{ii},tbin_edges_enc ); % get spike counts for each bin
    spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.
    spk_count = smoothdata(spk_count, 25, 'gaussian');
    
    data_h(:,ii) = spk_count;
end


% same for last 30mins sleep


    tstart = out.Sleep.csc.tvec(end)-60*60;
    tend =  out.Sleep.csc.tvec(end)-15*60; 

tbin_edges = tstart:cfg.bin_s:tend; % vector of time bin edges (for histogram)
tbin_centers = tbin_edges(1:end-1)+cfg.bin_s/2; % vector of time bin centers (for plotting)

data_r = [];
for ii = size(out.Sleep.S.t,2):-1:1
    
    
    spk_count = histc(out.Sleep.S.t{ii},tbin_edges); % get spike counts for each bin
    spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.
    spk_count = smoothdata(spk_count, 25, 'gaussian');
    
    data_r(:,ii) = spk_count;
end


[r, a] = ReactStrength(data_h, data_r);

[r_en] = ReactStrength(data_h, data_h);

%% plot the location of Aseemblies in encoding/recall
c_ord = winter(size(r_en,2))
figure(100)

subplot(1,2,1)
hold on
for aa = 1:size(r_en,2)

    [~,idx] = findpeaks(r_en(:,aa), 'minpeakheight', 5, 'minpeakdistance', 1/cfg.bin_s);

    idx_t = nearest_idx(tbin_centers_enc(idx), out.Encode.pos.tvec);
    plot(out.Encode.pos.data(1,idx_t), out.Encode.pos.data(2,idx_t), 'o', 'color', c_ord(aa,:))

end


    
   





%% try the sleep reactstrength


figure(101)
clf
ax(1) = subplot(5,1,1:3);
idx(1) = nearest_idx(tstart, out.Sleep.csc.tvec);
idx(2) = nearest_idx(tend, out.Sleep.csc.tvec);
plot(out.Sleep.csc.tvec(idx(1):idx(2)), out.Sleep.csc.data(1,idx(1):idx(2))*50000)

yyaxis right
plot(out.Sleep.csc.tvec(idx(1):idx(2)), out.Sleep.csc.data(2,idx(1):idx(2))*50000)

ax(2) = subplot(5,1,4:5);
plot(tbin_centers, r)

linkaxes(ax, 'x')


%% SWR times

swr = MS_SWR_detector(out.Sleep.csc, 'CSC21', 0); 
swr_cent = IVcenters(swr); 

react_t_avg = cell(size(r,2),1); 

win = .25; 
win = win/cfg.bin_s; 


for ii = length(swr_cent):-1:1
    idx = nearest_idx(swr_cent(ii), tbin_centers);

    for iR = size(r,2):-1:1

        react_t_avg{iR}(ii,:) = r((idx - win):(idx+win), iR);
    end
end


%% plot the SWR average for the assemblies. 

figure(102)
clf

hold on
 for iR = size(r,2):-1:1 
plot(-.25:cfg.bin_s:.25, mean(react_t_avg{iR}))

 end

 vline(0)
