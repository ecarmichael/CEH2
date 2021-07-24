function [wave_prop] = MS_get_wave_properties(S, wave,tvec,  plot_flag)
%% MS_get_wave_properties:
%
%
%
%    Inputs:
%    - S [struct]   Spike 'S' structure from LoadSpikes
%
%    - wave  [5 x nSamples]  waveform array.
%
%    - plot_flag: [logical]  0 = no plots, 1 = plots.
%
%    Outputs:
%    - wave_prop  [struct] containing:
%                    - burting_idx
%                    - spike_width
%                    - peak_val
%                    - trough_val
%                    - peak2trough
%
%
%
%
% EC 2021-07-22   initial version
%
%
%
%% initialize
if nargin < 4
    plot_flag = 0;
end

%%

spike_rate = length(S.t{1})./(tvec(end) - tvec(1));

log_isi = log10(mean(diff(S.t{1})));

%% get the autocorrelation (from vandermeerlab.org wiki: https://rcweb.dartmouth.edu/~mvdm/wiki/doku.php?id=analysis:nsb2016:week9)
xbin_centers = -.025-0.001:0.001:0.025+0.001; % first and last bins are to be deleted later
ac = zeros(size(xbin_centers));

for iSpk = 1:length(S.t{1})
    
    relative_spk_t = S.t{1} - S.t{1}(iSpk);
    
    ac = ac + hist(relative_spk_t,xbin_centers); % note that hist() puts all spikes outside the bin centers in the first and last bins! delete later.
    
end

xbin = xbin_centers(2:end-1); % remove unwanted bins
zero_idx = find(xbin == 0);
ac = ac(2:end-1);
ac(zero_idx) = 0;

% ac = ac./max(ac); % normalize

figure(900)
bar(xbin*1000, ac)
xlabel('Lag (ms)');


bin_idx_2_6 = zero_idx +2: zero_idx+6;
bin_idx_2_20 = zero_idx +2: zero_idx+20;

%% get spike waveform measures

wave_tvec = wave(1,:);

[~, best_wave_idx] = max((mean(abs(wave(2:end,:)),2)));

best_wave = wave(best_wave_idx+1, :);
[p_max, p_idx] = max(best_wave);
[t_max, t_idx] = min(best_wave);

if wave_tvec(t_idx)< wave_tvec(p_idx)
    best_wave = -best_wave;
    [p_max, p_idx] = max(best_wave);
    [t_max, t_idx] = min(best_wave);
end

pt_ratio = abs(p_max/t_max);
width = (wave_tvec(t_idx) - wave_tvec(p_idx));
if plot_flag
    figure(901)
    plot(wave_tvec,best_wave)
    
    hold on
    plot(wave_tvec(p_idx),p_max,"*")
    plot(wave_tvec(t_idx),t_max,"*")
    
    plot(wave_tvec([p_idx,t_idx]),[p_max,p_max],"-r")
    y_lim = ylim;
    ylim([y_lim(1) y_lim(2)*1.1]);
    plot(wave_tvec([t_idx,t_idx]),[p_max,t_max],"-r")
    text((wave_tvec(t_idx)-wave_tvec(p_idx))+wave_tvec(p_idx),p_max*1.1,num2str(width))
    
end
% width_arrow = annotation('doublearrow');
% width_arrow.Parent = width_arrow.CurrentAxes;  % associate annotation with current axes
% % now you can use data units
% width_arrow.X = [p_idx p_max];
% width_arrow.Y = [t_idx p_max];

%hold on
%x=[p_idx,p_max];
%y=[t_idx,p_max];
%annotation('doublearrow',x,y,'String','Width')
%hold off


wave_prop.burstingidx = [];
wave_prop.spike_width = width;
wave_prop.peak_val = p_max;
wave_prop.trough_val = t_max;
wave_prop.pt_ratio = pt_ratio;
