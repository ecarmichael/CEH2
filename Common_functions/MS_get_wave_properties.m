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

ISI = diff(S.t{1});


% ac = ac./max(ac); % normalize
if plot_flag
    figure(900)
    clf;
    subplot(2,3,3)
    bar(xbin*1000, ac)
    xlabel('Lag (ms)');
    
    subplot(2,3,6)
    histogram(log10(ISI),50,'Normalization', 'probability')
%     set(gca, 'XScale', 'log')
    set(gca, 'XTickLabel', 10.^get(gca, 'xtick'))
    xlabel('ISI (sec)');
    xline(prctile(log10(ISI), 10), '--', '10th prctl')
    xlim(log10([0.001 100]))
end

bin_idx_2_6 = zero_idx +2: zero_idx+6;
bin_idx_2_20 = zero_idx +2: zero_idx+20;

%% get spike waveform measures

wave_tvec = wave(1,:);

[~, best_wave_idx] = max((mean(abs(wave(2:end,:)),2)));

best_wave = wave(best_wave_idx+1, :);
[p_max, p_idx] = max(best_wave);
[t_max, ~] = min(best_wave(p_idx:end));
t_idx = find(best_wave == t_max); 

% % check that valley flows

% if wave_tvec(t_idx)< wave_tvec(p_idx)
%     best_wave = -best_wave;
%     [p_max, p_idx] = max(best_wave);
%     [t_max, t_idx] = min(best_wave);
% end


pt_ratio = abs(p_max/t_max);
width = (wave_tvec(t_idx) - wave_tvec(p_idx));

% properties from:
% get the valley.
[val_val, val_idx]= min(best_wave(p_idx:end)); % find the lowest point after the peak.

val_idx = val_idx + p_idx -1; % adjust for peak index.


% get the rise using the first positive slope ( T1 in
% http://humanspatialcognitionlab.org/wp-content/uploads/2014/10/ViskEtal07.pdf)

rise_idx = find(diff(best_wave)> 0, 1);

% fall idx using first negative after peak (T3)

fall_idx = find(best_wave(p_idx:end) <=0, 1);
fall_idx = fall_idx +p_idx-1; % compensate for max idx

% get the return to baseline ('0' or T5 in the Viskontas et al. 2007
% example).
return_idx = find(best_wave(val_idx:end) >=0, 1);
if isempty(return_idx)
    return_idx = length(best_wave(val_idx:end)); % if the wave doesnt return just get the last value
end
return_idx = return_idx+val_idx-1; %compensate for val_idx offset.


% get 'waveform duration' (ms);
wave_dur = wave_tvec(return_idx) - wave_tvec(rise_idx);

% get the slope ratio
slopes_ratio = abs((((best_wave(p_idx)) - best_wave(rise_idx)) / (wave_tvec(p_idx) - wave_tvec(rise_idx)))/...
    (((best_wave(fall_idx)) - best_wave(p_idx)) / (wave_tvec(fall_idx) - wave_tvec(p_idx))));

% get the 25% to 75% interval
rise_fall_inter =abs( ((wave_tvec(p_idx) - wave_tvec(rise_idx))*.25 +wave_tvec(rise_idx)) - ((wave_tvec(fall_idx) - wave_tvec(p_idx))*.75+wave_tvec(p_idx)));

% get the peak to valley ratio
pt_ratio = best_wave(p_idx) / best_wave(val_idx);

% bursting propensity
burst_idx = sum(ISI < .010) / sum(ISI > .010);
      
%%
      
if plot_flag
    figure(900)
    subplot(2,3,[1 2 4 5])
    plot(wave_tvec,best_wave)
    
    hold on



%         plot(wave_tvec([rise_idx,return_idx]),[t_max,t_max]-abs(t_max*.2),"-b")
        

    plot(wave_tvec([p_idx,t_idx]),[p_max,p_max],"-r")
    text(wave_tvec(t_idx), p_max - abs(t_max), 'Peak:Vally')
    y_lim = ylim;
    ylim([y_lim(1) y_lim(2)*1.1]);
    plot(wave_tvec([t_idx,t_idx]),[p_max,t_max],"-r")
    text((wave_tvec(t_idx)-wave_tvec(p_idx))+wave_tvec(p_idx),p_max*1.1,num2str(width))
    
    
        plot(wave_tvec(p_idx),p_max,"*")
    text(wave_tvec(p_idx),p_max,"peak")
    plot(wave_tvec(t_idx),t_max,"*")
    text(wave_tvec(t_idx),t_max,"valley")
    plot(wave_tvec(fall_idx),best_wave(fall_idx),"*")
    text(wave_tvec(fall_idx),best_wave(fall_idx), 'fall')
    plot(wave_tvec(return_idx),best_wave(return_idx),"*")
    text(wave_tvec(return_idx),best_wave(return_idx),"return")
    
    plot(wave_tvec(rise_idx),best_wave(rise_idx),"*")
    text(wave_tvec(rise_idx),best_wave(rise_idx),"rise")
    
    set(gca, 'XTickLabel', get(gca, 'xtick')*1000)
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
wave_prop.wave = wave;
wave_prop.firing_rate = spike_rate;
wave_prop.ISI = diff(S.t{1}); 
wave_prop.burstingidx = [];
wave_prop.spike_width = width;
wave_prop.peak_val = p_max;
wave_prop.trough_val = t_max;
wave_prop.pt_ratio = pt_ratio;
% Viskontas et al. 2006
wave_prop.wave_dur = wave_dur;
wave_prop.slopes_ratio = slopes_ratio;
wave_prop.rise_fall_inter = rise_fall_inter;
wave_prop.burst_idx = burst_idx;


wave_prop.auto_corr.ac = ac; 
wave_prop.auto_corr.xbin = xbin; 
