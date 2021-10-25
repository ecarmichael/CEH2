function tc = MS_speed_HD(cfg_in, S, pos,spd)
%% MS_speed_HD:
%
%
%
%    Inputs: 
%    -
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2021-10-16   initial version 
%
%
%
%% initialize

cfg_def = [];
cfg_def.spd_bins = .15;
cfg_def.hd_bins = 20;
cfg_def.interp = 1; % use interpolation
cfg_def.jump_thresh = 20; % how big a jump in the HD do you toleratore before removal. 
cfg = ProcessConfig2(cfg_def, cfg_in);

%%         
if cfg.interp
    keep_idx = (pos.data(3,:) ~= 0) & (pos.data(3,:) <= 360); % remove idx = 0.  Seems to be default state for NLX when unsure of dir.
    jumps = [0 rad2deg(circ_dist(deg2rad(pos.data(3,2:end)), deg2rad(pos.data(3,1:end-1))))]; % catch big jumps in HD angle
    keep_idx = keep_idx & abs(jumps) < cfg.jump_thresh;
    ulon=unwrap(pos.data(3,keep_idx)*pi/180)*180/pi;
    pos.data(3,:)= mod(interp1(pos.tvec(keep_idx),ulon,pos.tvec), 360);
end
%% get the hd by spd 
% flag bins that are outside of 0-360deg
keep_idx = (pos.data(3,:) ~= 0) & (pos.data(3,:) < 360);
keep_idx = keep_idx & spd.data > .2; 
% interp spikes to match speed bins

hd_bins = 0:cfg.hd_bins:360 - cfg.hd_bins;

spd_bins = .5:cfg.spd_bins:floor(max(spd.data)/.8);
vec_occ = nan(length(hd_bins)-1,length(spd_bins)-1);
spk_vec = nan(size(vec_occ));



spk_ang = interp1(pos.tvec,pos.data(3,:),S.t{1},'linear');


for ihd = 1:length(hd_bins)-1

    if ihd == 1
        idx = (pos.data(3,:) <=hd_bins(ihd)+cfg.hd_bins);
    elseif ihd == length(hd_bins)-1
        idx = (pos.data(3,:) > hd_bins(ihd)+cfg.hd_bins);
    else
        idx =(hd_bins(ihd) < pos.data(3,:)) & (pos.data(3,:) <= hd_bins(ihd)+cfg.hd_bins);
    end
    if sum(idx & keep_idx) ~= 0
        vec_occ(ihd,:) = histcounts(spd.data(idx & keep_idx), spd_bins);
        
        
        spk_int = interp1(spd.tvec(idx & keep_idx),spd.data(idx & keep_idx),S.t{1},'linear');
        
        spk_vec(ihd,:) = histcounts(spk_int, spd_bins);
    end

end

vec_occ(vec_occ < 10) = NaN;

tc = spk_vec./vec_occ.* mode(diff(spd.tvec));

tc(isinf(tc)|isnan(tc)) = 0; 

% smooth with gaussian if you like

if cfg.smooth
   vec_occ = imgaussfilt(vec_occ,cfg.smooth_sd) ;
   spk_vec = imgaussfilt(spk_vec,cfg.smooth_sd) ;
   tc = imgaussfilt(tc,cfg.smooth_sd) ;

    
end

%% compute the mean vector length vs shuffle

occ_h = histcounts(pos.data(3,:),360/cfg.hd_bins);
spk_h = histcounts(spk_ang,360/cfg.hd_bins);
tc_h = spk_h./ occ_h .* mode(diff(spd.tvec));


circ_r(tc_h)

%% test figure

figure
subplot(2,3,1)
title('Occ')
polarplot(deg2rad(hd_bins),occ_h);


subplot(2,3,2)
title('spk')
polarplot(deg2rad(hd_bins),spk_h);

subplot(2,3,3)
title('tc')

polarplot(deg2rad(hd_bins),(tc_h));


subplot(2,3,4)
imagesc([hd_bins 360],spd_bins, vec_occ');
set(gca, 'xtick', [0 180 360], 'ytick', 0:floor(max(spd.data)/.8));
ylabel('velocity (cm/s)'); xlabel('HD (deg)');

axis xy
colormap([0,0,0; parula])
caxis([0 max(vec_occ, [], 'all')])

subplot(2,3,5)
title('Spk')
imagesc([hd_bins 360],spd_bins, spk_vec');
set(gca, 'xtick', [0 180 360], 'ytick', 0:floor(max(spd.data)/.8));
ylabel('velocity (cm/s)'); xlabel('HD (deg)');
axis xy
colormap([0,0,0; parula])
caxis([0 max(spk_vec, [], 'all')])


subplot(2,3,6)
title('tc')
imagesc([hd_bins 360],spd_bins, tc');
set(gca, 'xtick', [0 180 360], 'ytick', 0:floor(max(spd.data)/.8))
ylabel('velocity (cm/s)'); xlabel('HD (deg)');
axis xy
colormap([0,0,0; parula])
caxis([0 max(tc, [], 'all')/2])

