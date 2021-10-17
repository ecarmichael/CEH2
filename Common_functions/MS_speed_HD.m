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
cfg_def.spd_bins = .25;
cfg_def.hd_bins = 10;

cfg = ProcessConfig2(cfg_def, cfg_in);



%% get the hd by spd 
% flag bins that are outside of 0-360deg
keep_idx = (pos.data(3,:) ~= 0) & (pos.data(3,:) < 360);
keep_idx = keep_idx & spd.data > .2; 
% interp spikes to match speed bins


spd_bins = 0:cfg.spd_bins:5;
vec_occ = nan(cfg.hd_bins-2,length(spd_bins)-1);
spk_vec = nan(size(vec_occ));


hd_bins = 0:cfg.hd_bins:360 - cfg.hd_bins;

for ihd = 1:length(hd_bins)

    if ihd == 1
        idx = (pos.data(3,:) <=hd_bins(ihd)+cfg.hd_bins);
    elseif ihd == length(hd_bins)
        idx = (pos.data(3,:) > hd_bins(ihd)+cfg.hd_bins);
    else
        idx =(hd_bins(ihd) < pos.data(3,:)) & (pos.data(3,:) <= hd_bins(ihd)+cfg.hd_bins);
    end
    if sum(idx & keep_idx) ~= 0 
    vec_occ(ihd,:) = histcounts(spd.data(idx & keep_idx), spd_bins);
    
    vec_occ = vec_occ .* mode(diff(spd.tvec)); 

    spk_int = interp1(spd.tvec(idx & keep_idx),spd.data(idx & keep_idx),S.t{1},'linear'); 

    spk_vec(ihd,:) = histcounts(spk_int, spd_bins);
    end

end

tc = vec_occ./spk_vec;

tc(isinf(tc)|isnan(tc)) = 0; 
%% test figure

figure
imagesc(hd_bins, spd_bins, tc');
set(gca, 'xtick', [0 180 360], 'ytick', [0:5])
axis xy
colormap([0,0,0; parula])
caxis([0 max(tc, [], 'all')])




