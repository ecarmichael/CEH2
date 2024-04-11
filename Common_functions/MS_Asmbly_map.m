function [map_out, place_idx] = MS_Asmbly_map(P_cells, place, min_N_place)
%% MS_Asmbly_map:
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
% EC 2024-01-24   initial version 
%
%
%
%% initialize


%% loop over assemblies and interpolate the spatial maps of the cells

for ii = length(P_cells):-1:1

      % create a place map for the assembly;
      P_map{ii} = [];  P_cent{ii} = []; no_map_idx = [];
      
%       fprintf('Assembly # %.0f had %.0f place cells\n', ii, sum(place.is(P_cells{ii})))
      
      for jj = 1:length(P_cells{ii})
          
          place_int = interp1(place.p_bins,place.map(P_cells{ii}(jj),:),  place.p_bins(1):1:place.p_bins(end));
          
          if sum(place_int) == 0; no_map_idx(jj) = 0; else no_map_idx(jj) =1; end
          P_map{ii}(jj,:) =  place_int;
          P_cent{ii}(jj) = place.centroids(P_cells{ii}(jj));
      end

      
      for iS = 500:-1:1
          shuff_var(iS) = var(randsample(place.centroids, length(P_cent{ii})));
          
          % with the max position
          [~, these_peaks] = max(place.map(randsample(find(place.is),length(no_map_idx)),:),[],2);
          
          shuff_peak(iS) = var(these_peaks);
      end
      
      [~, P_peak{ii}] = max(place.map(P_cells{ii},:), [],2);
      
      P_peak_z(ii) = (var(P_peak{ii}) - mean(shuff_peak))/std(shuff_peak);
      
      P_cent_z(ii) = (var(P_cent{ii}) - mean(shuff_var))/std(shuff_var);
      
      map_out{ii} = []; 
      map_out{ii}.bins = place.p_bins(1):1:place.p_bins(end);
      map_out{ii}.map = P_map{ii}; 
      map_out{ii}.map_mean = nanmean(P_map{ii},1);
      map_out{ii}.cent = P_cent{ii}'; 
      map_out{ii}.peak = P_peak{ii};
      map_out{ii}.peak_z = P_peak_z(ii); 
      map_out{ii}.cent_z = P_cent_z(ii); 
      map_out{ii}.place_idx = place.is(P_cells{ii}); 
      map_out{ii}.MI = place.MI(P_cells{ii});
      map_out{ii}.rate = place.peak_rate(P_cells{ii}); 
      map_out{ii}.bin_size = place.bin_s; 

      place_idx(ii) = logical(sum(map_out{ii}.place_idx) >= min_N_place);

end



