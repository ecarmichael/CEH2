%% sandbox_PCA_ICA_CA

code_dir = 'C:\Users\ecarm\Downloads\Dos-Santos Assembly ICA\Dos-Santos Assembly ICA';

RnR_dir = 'C:\Users\ecarm\Documents\GitHub\RnR_methods';

p_cell_dir = 'C:\Users\ecarm\Documents\GitHub\Replay_REM\space_suite-pval';



addpath(p_cell_dir);
addpath(genpath(RnR_dir));
addpath(code_dir);

%%  load some data

load('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Radial\Rad_maze\dSub_g8_E\customEntValHere\2023_06_02\ms.mat')


e = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Radial\Rad_maze\dSub_g8_E\customEntValHere\2023_06_02\12_36_21\My_WebCam';
r = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Radial\Rad_maze\dSub_g8_E\customEntValHere\2023_06_02\17_19_09\My_WebCam';

enc = [];
rec = [];
cfg_pos = [];
cfg_pos.conv_fac = [5.98 5.65]; 
[~, enc_behav] = MS_DLC2TSD(e,[], cfg_pos.conv_fac); 
[~, rec_behav] = MS_DLC2TSD(r,[],  cfg_pos.conv_fac); 


%% pull out the encoding and recall phases

enc_idx = [1 length(ms.tvecs{1})]; 
enc.time = ms.tvecs{1}; 

enc.Binary = ms.Binary(enc_idx(1):enc_idx(2),:); 
enc.RawTraces = ms.RawTraces(enc_idx(1):enc_idx(2),:); 
enc.deconv = ms.deconv(enc_idx(1):enc_idx(2),:); 
enc.denoise = ms.denoise(enc_idx(1):enc_idx(2),:); 


rec_idx = [length(ms.time)+1 - length(ms.tvecs{end}) length(ms.time)]; 
rec.time = ms.tvecs{end}; 
rec.RawTraces = ms.RawTraces(rec_idx(1):rec_idx(2),:); 
rec.Binary = ms.Binary(rec_idx(1):rec_idx(2),:); 
rec.deconv = ms.deconv(rec_idx(1):rec_idx(2),:); 
rec.denoise = ms.denoise(rec_idx(1):rec_idx(2),:); 

% algin behaviour
enc_behav = MS_align_data(enc_behav, enc);
rec_behav = MS_align_data(rec_behav, rec);

enc_move_idx = enc_behav.speed > 3; 
rec_move_idx = rec_behav.speed > 3; 

%% extract place cell metrics
%     PARAMS.data.ca_time= enc.time;
%     PARAMS.data.ca_data=enc.RawTraces ;
%     PARAMS.data.behav_time=enc_behav.time;
%     PARAMS.data.behav_vec=enc_behav.position(:,1);
%     PARAMS.data.num_surrogates=1000;
%     
% [PC_properties]= Bayesian_JS2(inter_dir, info, PARAMS, ms, behav)
%% try the assembly script

Csp = enc.deconv./enc.denoise;
Csp = Csp > 0.01;
enc.Csp = Csp;

Csp = rec.deconv./rec.denoise;
Csp = Csp > 0.01;
rec.Csp = Csp;

%% try again with histc
binsize = .5;
tbin_edges = enc.time(1):binsize:enc.time(end); % vector of time bin edges (for histogram)
tbin_centers = tbin_edges(1:end-1)+binsize/2; % vector of time bin centers (for plotting)

data_e = [];
for ii = size(Csp,2):-1:1
    this_cell = enc.time(find(enc.Csp(: ,ii)));
    spk_count = histc(this_cell,tbin_edges); % get spike counts for each bin
    spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.
    
    data_e(:,ii) = spk_count;
end

tvec_e = tbin_centers;


%% try again with histc
binsize = .5;
tbin_edges = rec.time(1):binsize:rec.time(end); % vector of time bin edges (for histogram)
tbin_centers_r = tbin_edges(1:end-1)+binsize/2; % vector of time bin centers (for plotting)

data_r = [];
for ii = size(Csp,2):-1:1
    this_cell = rec.time(find(rec.Csp(: ,ii) ));
    spk_count = histc(this_cell,tbin_edges); % get spike counts for each bin
    spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.
    
    data_r(:,ii) = spk_count;
end

tvec_r = tbin_centers_r;

%% try the assembly code....
rng default
Rec_Temp = assembly_patterns(data_e');

rng default
time_proj_e = assembly_activity(Rec_Temp,data_e');
rng default
time_proj_r = assembly_activity(Rec_Temp,data_r');



%% plot check
c_ord = parula(50);
n = 46;
figure(101)
clf
ax(1) = subplot(10,1,1);
plot(enc_behav.time, enc_behav.position(:,1:2)); 
yyaxis right
plot(enc_behav.time, enc_behav.speed); 
xlim([enc.time(1) enc.time(end)])

ax(2) = subplot(10,1,2:4); 
hold on
for ii = 1:n
plot(enc.time, zscore(enc.RawTraces(:,ii))+ii*10, 'color', c_ord(ii,:)); 
plot(enc.time, enc.Binary(:,ii)*4+ii*10, 'color', c_ord(ii,:)); 
end
xlim([enc.time(1) enc.time(end)])

ax(3) = subplot(10, 1, 5);
plot(tvec_e, time_proj_e)
linkaxes(ax, 'x')
xlim([enc.time(1) enc.time(end)])

ar(1) = subplot(10,1,6);
plot(rec_behav.time, rec_behav.position(:,1:2)); 
yyaxis right
plot(rec_behav.time, rec_behav.speed); 

ar(2) = subplot(10,1,7:9);
hold on
for ii = 1:n
plot(rec.time, zscore(rec.RawTraces(:,ii))+ii*10, 'color', c_ord(ii,:)); 
plot(rec.time, rec.Binary(:,ii)*4+ii*10, 'color', c_ord(ii,:)); 
end
ar(3) = subplot(10, 1, 10);
plot(tvec_r, time_proj_r)

linkaxes(ar, 'x')
xlim([rec.time(1) rec.time(end)])

%% keep assemblies with strong activations
keep_idx = zeros(1,size(Rec_Temp,2));
Ass_pos_cells = cell(size(keep_idx));
for ii = size(Rec_Temp,2):-1:1
    
    if max(Rec_Temp(:, ii)) > 0.2
        keep_idx(ii) = 1;
        
        z_weight = zscore(Rec_Temp(:,ii));
        
        Ass_pos_cells{ii} = find(z_weight > 2);
        
    end
    
    
end

Ass_pos = Rec_Temp;
Ass_pos(:,~keep_idx) = [];
Ass_pos_cells(~keep_idx) = [];
Ass_idx = find(keep_idx);

time_proj_pos = time_proj_r;
time_proj_pos(~keep_idx,:) = [];

% Ass_z = (size(Ass_pos,2) - mean(Ass_shuff))/std(Ass_shuff);
% fprintf('%.0f Positive Assemblies detected. Chance level is %.1f. zscore = %.1fSD\n', size(Ass_pos,2), mean(Ass_shuff), Ass_z)

%%
%% stem plot ensembles
figure(301);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
figure(302);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

Ass_map = cell(size(Ass_pos,2),1);
Ass_mean_map = [];
Ass_pcells = Ass_map;
c_ord = parula(size(Ass_pos,2)+2);

cmap = parula(256);
cmap(1,:) = 0;
for ii = size(Ass_pos,2):-1:1
    figure(301)
    subplot(4, ceil(size(Ass_pos,2)/4),ii)
    hold on
    stem(Ass_pos(:,ii), 'color', c_ord(ii,:))
    view(90,90)
    
    stem(Ass_pos_cells{ii}, Ass_pos(Ass_pos_cells{ii},ii), 'color', c_ord(ii,:), 'MarkerFaceColor', c_ord(ii,:))
    title(['Assembly #' num2str(ii)])
    
%     
%     figure(302)
%     subplot(4, ceil(size(Ass_pos,2)/4),ii)
%     hold on
%     this_ass_map = []; these_place = zeros(size(Ass_pos_cells{ii}));
%     fprintf('Assembly # %.0f had %.0f place cells\n', ii, sum(place.is(Ass_pos_cells{ii})))
%     for jj = 1:length(Ass_pos_cells{ii})
%         
%         if place.is(Ass_pos_cells{ii}(jj))
%             these_place(jj) = 1;
%             place_int = interp1(p_bins,place.map(Ass_pos_cells{ii}(jj),:),  p_bins(1):1:p_bins(end));
%             
%             this_ass_map = [this_ass_map ; place_int];
%         end
%     end
%     
%     if ~isempty(this_ass_map)
%         imagesc(p_bins(1):1:p_bins(end), 1:length(sum(these_place)),  this_ass_map)
%         set(gca,'ytick',1:size(this_ass_map,1), 'YTickLabel',  Ass_pos_cells{ii}(logical(these_place))');
%         xlim([p_bins(1) p_bins(end)])
%         ylim([.5 size(this_ass_map,1)+.5])
%         colormap(cmap)
%         Ass_map{ii} = this_ass_map;
%         Ass_pcells{ii} = Ass_pos_cells{ii}(logical(these_place));
%         Ass_mean_map(ii,:) = mean(this_ass_map,1);
%     else
%         
%         Ass_map{ii} = NaN(size(p_bins(1):1:p_bins(end)));
%         Ass_pcells{ii} = NaN;
%         Ass_mean_map(ii,:) = NaN(size(p_bins(1):1:p_bins(end)));
%     end
%     title(['Assembly #' num2str(ii)])
    
    
    
    %     ms_t = MS_append_sharp_SFPs(ms_trk);
    
    
    %     MS_plot_all_SFPs(imgaussfilt(ms_t.SFPs_sharp(:,:, Ass_pos_cells{ii}),2))
    
end

%% get the assembly triggered position average.
win = floor(2.5 * mode(diff(enc_behav.time)));

x_bins = 10:2.5:100;
y_bins = 0:2.5:80; 
kernel = gausskernel([1 1],2); % 2d gaussian in bins

% compute occupancy encode
occ_hist = histcn(rec_behav.position(rec_move_idx,1:2),x_bins,y_bins); % 2-D version of histc()
occ_hist = conv2(occ_hist,kernel,'same');

occ_hist = occ_hist .* (1/30); % convert samples to seconds using video frame rate (30 Hz)
% no_occ_idx = find(occ_hist < 5); % NaN out bins never visited
% occ_hist(no_occ_idx) = NaN;

figure(5)
clf
n = 5;%ceil(size(time_proj_r,1)/2);
m = 2;

nFig = 5; count = 0; 
for ii = 1: size(time_proj_r,1)
        count = count+1; 

    if count > nFig
        figure(ii)
        count = 1; 
    end
    subplot(m,n,count)
    hold on
    [~, p_idx] = findpeaks(zscore(time_proj_r(ii,:)),'MinPeakHeight', 1.96 ,'MinPeakDistance', 2*floor(mode(diff(rec.time))));
    
    this_pos = [];
    plot(rec_behav.position(:, 1), rec_behav.position(:,2), '.k')
    for ip = 1:length(p_idx)
        this_idx = nearest_idx(tbin_centers_r(p_idx(ip)), rec_behav.time);
        plot(rec_behav.position(this_idx, 1), rec_behav.position(this_idx,2), 'or')
%         pos_mat(:,:,ip) = 
%         if ((this_idx - win) >=0) && ((this_idx + win)<= length(enc_behav.time))
%             this_pos(ip,:) = enc_behav.position(this_idx - win:this_idx+win,1);
%             plot((-win:win)/mode(diff(behav.time)), this_pos(ip,:), 'color',[c_ord(ii+4,:) .5])
%         end
    end
    spk_x = interp1(rec_behav.time, rec_behav.position(:,1),tbin_centers_r(p_idx)','linear');
    spk_y = interp1(rec_behav.time, rec_behav.position(:,2),tbin_centers_r(p_idx)','linear');
    
     % compute the rate map
    spk_hist = histcn([spk_x, spk_y],x_bins,y_bins);
    spk_hist = conv2(spk_hist,kernel, 'same');
    tc{ii} = spk_hist./occ_hist;
    subplot(m,n, count+n)
    imagesc(x_bins, y_bins, tc{ii})
    set(gca, 'YDir', 'normal')
    
    
%     plot((-win:win)/mode(diff(behav.time)), mean(this_pos), 'color',[c_ord(ii+4,:) 1], 'linewidth', 3)
%     xlim([-win/mode(diff(behav.time)) win/mode(diff(behav.time))]);
%     %         set(gca, 'color', 'k')
    title(['Assembly #' num2str(ii)])

end

%% Colorcode the assemblies 
figure(401);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);





