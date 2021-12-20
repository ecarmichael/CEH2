%% waveform classification sandbox

% add paths to code bases here:




%% load some data from the inter mediate folder. 

inter_dir = '/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/dSubiculum/inter/All_maze_cells';

% get a list of cells;

inter_cells = dir('M*'); 

for iCell = 1%:length(inter_cells) % loop over all available cells. 
    
    fprintf('Loading cell: %s...\n', inter_cells(iCell).name); % print the name of the session being loaded (optional)
    
    load(inter_cells(iCell).name); % laod the intermediate file containing the spiking and waveform data. 
    
    
    
    %% get the wave form properties here
    % get the waves across channels.  The first column is time, then electrode channels. 
    wave_t = This_Cell.wave.wave(:,1); 
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
       

   if  max_idx < neg_max_idx
       % if the negative max_index is less than the normal max idx (as
       % would be the case with an inverted spike, then use the negative
       % values. 
       fprintf('Spike appears to be inverted, inverting back...\n')
       waves = -waves;
       [max_val, max_idx] = max(waves(:, best_chan));
   end
   
      % get the valley. 
      [~, val_idx]= min(waves(max_idx:end, best_chan)); % find the lowest point after the peak. 

      val_idx = val_idx + max_idx -1; % adjust for peak index. 
      
      % plot everything to double check.
      figure(101)
      clf
      plot(wave_t, waves, 'linewidth', 2);
      hold on
      plot(wave_t(max_idx), waves(max_idx, best_chan), '+r', 'markersize',12,'linewidth', 2);
      plot(wave_t(val_idx), waves(val_idx, best_chan), 'xb', 'markersize',12, 'linewidth', 2);
      legend({'Chan 1', 'Chan 2', 'Chan 3', 'Chan 4','peak', 'trough'});
      
   %% compute the time differences between the peak and the valley.  
   
   % get the peak to valley ratio
   pvr = max_val / waves(val_idx, best_chan)
   
   % get the peak to valley time
   
   
   % get symetry
   
   
   
   
end


%% try some population level plotting and maybe a classification method. 



for iCell = 1:length(inter_cells) % loop over all available cells. 
    
    fprintf('Loading cell: %s...\n', inter_cells(iCell).name); % print the name of the session being loaded (optional)
    
    load(inter_cells(iCell).name); % laod the intermediate file containing the spiking and waveform data. 
    
start_t = This_Cell.csc.tvec(find(diff(This_Cell.csc.tvec) > 1)+2); 
end_t = This_Cell.csc.tvec(find(diff(This_Cell.csc.tvec) > 1)); 


r_S = restrict(This_Cell.S, start_t(1),end_t(2)); 

    spike_stats.p_t(iCell) = This_Cell.wave.pt_ratio; 
    spike_stats.ISI(iCell) = mode(diff(r_S.t{1}))*1000;
    spike_stats.rate(iCell) = length(r_S.t{1}) / (end_t(2) - start_t(1)); 
    spike_stats.peak_val(iCell) = This_Cell.wave.peak_val; 
    spike_stats.spike(iCell) = This_Cell.wave.spike_width; 

    
end

%% plot a 3d scatter of the spike properties


figure(1010)
scatter(spike_stats.rate, spike_stats.ISI)