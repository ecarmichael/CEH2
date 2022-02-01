%% waveform classification sandbox

% add paths to code bases here:




%% load some data from the inter mediate folder. 

    addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'));

% inter_dir = '/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/dSubiculum/inter/All_maze_cells';
inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\dSubiculum\inter\All_maze_cells';

% get a list of cells;

inter_cells = dir('M*'); 

%%


for iCell = length(inter_cells):-1:1 % loop over all available cells. 
    
    fprintf('Loading cell: %s...\n', inter_cells(iCell).name); % print the name of the session being loaded (optional)
    
    load(inter_cells(iCell).name); % laod the intermediate file containing the spiking and waveform data. 
    
    
    
    %% get the wave form properties here
    % get the waves across channels.  The first column is time, then electrode channels. 
    wave_t = 0:0.06:1.94; % NLX spike time window in 32 samples. 
    waves = This_Cell.wave.wave(:,2:5); 
    
    for ii =1:4
        if sum(waves(:, ii) ==0) == length(waves(:,ii)); % check if this column is all zeros. 
            waves(:,ii) = NaN; % if it is, replace it with NaN ('Not a number') which will not show up in plots or stats. 
        end
    end % end the loop over channels. 
    
%     % plot to have a look
%     figure(100)
%     plot(wave_t, waves);
%     legend({'Chan 1', 'Chan 2', 'Chan 3', 'Chan 4'}); 
    
    %% get some wave properties
   [~, best_chan] = max(max(abs(waves),[],1)); % find the column index for the channel with the max value in either positive or negative directions. Use this for computing properties. 
    
   % check is the waveform is inverted. 
       [max_val, max_idx] = max(waves(:, best_chan)); 

       [neg_max_val, neg_max_idx] = max(-waves(:, best_chan)); 
       

   if  max_idx > neg_max_idx
       wave_dur(iCell) = NaN; 
       slopes_ratio(iCell) = NaN; 
       rise_fall_inter(iCell) = NaN; 
        pt_ratio(iCell) = NaN; 
        burst_idx(iCell) = NaN;
        firing_rate(iCell) = NaN; 
       continue
       % if the negative max_index is less than the normal max idx (as
       % would be the case with an inverted spike, then use the negative
       % values. 
%        fprintf('Spike appears to be inverted, inverting back...\n')
%        waves = -waves;
%        [max_val, max_idx] = max(waves(:, best_chan));
   end
   
      % get the valley. 
      [~, val_idx]= min(waves(max_idx:end, best_chan)); % find the lowest point after the peak. 

      val_idx = val_idx + max_idx -1; % adjust for peak index. 
      
      
      % get the rise using the first positive slope ( T1 in 
      % http://humanspatialcognitionlab.org/wp-content/uploads/2014/10/ViskEtal07.pdf)
      
      rise_idx = find(diff(waves(:,best_chan))> 0, 1); 
      
      % fall idx using first negative after peak (T3)
      
      fall_idx = find(waves(max_idx:end, best_chan) <=0, 1);
      fall_idx = fall_idx +max_idx-1; % compensate for max idx
      
      % get the return to baseline ('0' or T5 in the Viskontas et al. 2007
      % example). 
      return_idx = find(waves(val_idx:end,best_chan) >=0, 1);
      if isempty(return_idx)
          return_idx = length(waves(val_idx:end,best_chan)); % if the wave doesnt return just get the last value
      end
      return_idx = return_idx+val_idx-1; %compensate for val_idx offset. 
      
     
%       % plot everything to double check.
%       figure(101)
%       clf
%       plot(wave_t, waves, 'linewidth', 2);
%       hold on
%       hline(0)
%       plot(wave_t(max_idx), waves(max_idx, best_chan), '+r', 'markersize',12,'linewidth', 2);
%       plot(wave_t(val_idx), waves(val_idx, best_chan), 'xb', 'markersize',12, 'linewidth', 2);
%       plot(wave_t(rise_idx), waves(rise_idx, best_chan), 'sg', 'markersize',12, 'linewidth', 2);
%       plot(wave_t(fall_idx), waves(fall_idx, best_chan), 'om', 'markersize',12, 'linewidth', 2);
%       plot(wave_t(return_idx), waves(return_idx, best_chan), 'dc', 'markersize',12, 'linewidth', 2);
%  
%       legend({'Chan 1', 'Chan 2', 'Chan 3', 'Chan 4','peak', 'trough', 'rise', 'fall', 'return'});
%       
      
      %
   %% compute the time differences between the peak and the valley.  
   
   % get 'waveform duration' (ms); 
   wave_dur(iCell) = wave_t(return_idx) - wave_t(rise_idx); 
   
   % get the slope ratio
   slopes_ratio(iCell) = abs((((waves(max_idx,best_chan)) - waves(rise_idx,best_chan)) / (wave_t(max_idx) - wave_t(rise_idx)))/...
       (((waves(fall_idx,best_chan)) - waves(max_idx,best_chan)) / (wave_t(fall_idx) - wave_t(max_idx)))); 
   
   % get the 25% to 75% interval
   rise_fall_inter(iCell) =abs( ((wave_t(max_idx) - wave_t(rise_idx))*.25 +wave_t(rise_idx)) - ((wave_t(fall_idx) - wave_t(max_idx))*.75+wave_t(max_idx))); 
   
   % get the peak to valley ratio
   pt_ratio(iCell) = waves(max_idx,best_chan) / waves(val_idx, best_chan); 
   
   % bursting propensity
   ISI = diff(This_Cell.S.t{1}); 
   burst_idx(iCell) = sum(ISI < .010) / sum(ISI > .010);
   
   % get the firing rate corrected for start and stops. 
   start_t = This_Cell.csc.tvec(find(diff(This_Cell.csc.tvec) > 1)+2);
   end_t = This_Cell.csc.tvec(find(diff(This_Cell.csc.tvec) > 1));
   
   firing_rate(iCell) = length(This_Cell.S.t{1}) / ((end_t(1) - This_Cell.csc.tvec(1)) + (end_t(2) - start_t(1)) + (end_t(3) - start_t(2)) +(This_Cell.csc.tvec(end) - start_t(3))); 
   
   
end

%% plot some properties in 3D

% figure(101)
% scatter3(rise_fall_inter, burst_idx, firing_rate)
firing_rate = firing_rate(~isnan(firing_rate)); 
burst_idx = burst_idx(~isnan(burst_idx)); 
rise_fall_inter = rise_fall_inter(~isnan(rise_fall_inter)); 


% k-means
data_in = [firing_rate', burst_idx', rise_fall_inter'];

[g_idx, n_idx, frq_idx]= MS_kmean_scatter(data_in, 3,[3,2,1], 50);
clus1 = frq_idx(1); % get the top two clusters. 
clus2 = frq_idx(2); 

close(gcf);

figure(221)
subplot(2,3,[1,2,4,5])
hold on
col_ord = linspecer(2); 
scatter3(rise_fall_inter(g_idx==clus1), burst_idx(g_idx==clus1), firing_rate(g_idx==clus1), 50,col_ord(1,:), 'filled');
scatter3(rise_fall_inter(g_idx==clus2), burst_idx(g_idx==clus2), firing_rate(g_idx==clus2), 50,col_ord(2,:), 'filled');

grid on
xlabel('Rise-fall interval (ms)'); ylabel('Burst ISI ratio'); zlabel('Firing rate (Hz)')
legend('Cluster 1', 'Cluster 2')

subplot(2,3,3)
hold on
b = bar([1,2], [mean(firing_rate(g_idx==clus1)) mean(firing_rate(g_idx==clus2))]');
b.FaceColor = 'flat';
b.CData(1,:) =col_ord(1,:); 
b.CData(2,:) =col_ord(2,:); 

eb = errorbar([1,2], [mean(firing_rate(g_idx==clus1)) mean(firing_rate(g_idx==clus2))],...
    [mean(firing_rate(g_idx==clus1))-(std(firing_rate(g_idx==clus1))/sqrt(length(firing_rate(g_idx==clus1)))),...
     mean(firing_rate(g_idx==clus2))-(std(firing_rate(g_idx==clus2))/sqrt(length(firing_rate(g_idx==clus2))))],...
        [mean(firing_rate(g_idx==clus1))+(std(firing_rate(g_idx==clus1))/sqrt(length(firing_rate(g_idx==clus1)))),...
        mean(firing_rate(g_idx==clus2))+(std(firing_rate(g_idx==clus2))/sqrt(length(firing_rate(g_idx==clus2))))]);
    eb.Color = [0 0 0];
    
    ylabel('Firing rate (Hz)'); 
    
subplot(2,3,6)
hold on
b =bar([1,2], [mean(burst_idx(g_idx==clus1)) mean(burst_idx(g_idx==clus2))]);
b.FaceColor = 'flat';
b.CData(1,:) =col_ord(1,:); 
b.CData(2,:) =col_ord(2,:); 
eb = errorbar([1,2], [mean(burst_idx(g_idx==clus1)) mean(burst_idx(g_idx==clus2))],...
    [mean(burst_idx(g_idx==clus1))-(std(burst_idx(g_idx==clus1))/sqrt(length(burst_idx(g_idx==clus1)))),...
     mean(burst_idx(g_idx==clus2))-(std(burst_idx(g_idx==clus2))/sqrt(length(burst_idx(g_idx==clus2))))],...
        [mean(burst_idx(g_idx==clus1))+(std(burst_idx(g_idx==clus1))/sqrt(length(burst_idx(g_idx==clus1)))),...
        mean(burst_idx(g_idx==clus2))+(std(burst_idx(g_idx==clus2))/sqrt(length(burst_idx(g_idx==clus2))))]);
    
    eb.Color = [0 0 0];
    ylabel('Burst ratio')

%% try some population level plotting and maybe a classification method. 



% for iCell = 1:length(inter_cells) % loop over all available cells. 
%     
%     fprintf('Loading cell: %s...\n', inter_cells(iCell).name); % print the name of the session being loaded (optional)
%     
%     load(inter_cells(iCell).name); % laod the intermediate file containing the spiking and waveform data. 
%     
% start_t = This_Cell.csc.tvec(find(diff(This_Cell.csc.tvec) > 1)+2); 
% end_t = This_Cell.csc.tvec(find(diff(This_Cell.csc.tvec) > 1)); 
% 
% 
% r_S = restrict(This_Cell.S, start_t(1),end_t(2)); 
% 
%     spike_stats.p_t(iCell) = This_Cell.wave.pt_ratio; 
%     spike_stats.ISI(iCell) = mode(diff(r_S.t{1}))*1000;
%     spike_stats.rate(iCell) = length(r_S.t{1}) / (end_t(2) - start_t(1)); 
%     spike_stats.peak_val(iCell) = This_Cell.wave.peak_val; 
%     spike_stats.spike(iCell) = This_Cell.wave.spike_width; 
% 
%     
% end

% %% plot a 3d scatter of the spike properties
% 
% 
% figure(1010)
% scatter(spike_stats.rate, spike_stats.ISI)