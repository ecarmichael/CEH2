function overMap = MS_map_overlap(maps1, maps2,plot_flag,  names)
%% MS_map_overlap: Compares the spatial overlap between two sets of rate maps. Only includes cells that are present in both the maps1 and maps2.
%
%
%
%    Inputs:
%    - maps1: [struct]    output from MS_get_place_field
%
%    - maps2: [struct]    maps for comparison
%
%    - plot_flag: [logical]     1 for plots, 0 for skip.
%
%    - names: [2 x 1 cell array]  contains strings for the names of the
%    maps. exmaple {'Encoding', 'Recall'}
%
%    Outputs:
%    - overMap: [strcut]    Contains spatial overlap measures for each cell
%    pair between maps1 and maps2
%
%
%
%
% EC 2023-08-18   initial version
%
%
%
%% initialize

if nargin < 3
    plot_flag = 1;
    names = {'map 1', 'map 2'};
elseif nargin < 4
    names = {'map 1', 'map 2'};
    
end


%% get the pairs of cells for maps1 and maps 2
keep_idx = ~cellfun('isempty', maps1.tc) & ~cellfun('isempty', maps2.tc);

cell_ids = maps1.label(keep_idx);


%% loop over cells and get the 2d spatial correlation
overMap = [];

for iC = 1:length(cell_ids)
    
    overMap.label{iC} = cell_ids{iC};
    this_map1 = maps1.tc{find(ismember(maps1.label, cell_ids{iC}))};
    this_map1(isnan(this_map1)) = 0;
    
    this_map2 = maps2.tc{find(ismember(maps2.label, cell_ids{iC}))};
    this_map2(isnan(this_map2)) = 0;
    
    overMap.corr(iC) = corr2(this_map1, this_map2);
    
end



%% plot the output if asked


if plot_flag
    
    %subplot grid.
    m = 4;
    n = 4;
    s1_idx = 3:2:(n*m);
    s2_idx = 4:2:(n*m);
    
    % plot
    figure(202)
    clf
    ip = 0;
    
    subplot(m, n, 1)
    imagesc(maps1.occ);
    title(names{1})
    set(gca, 'xtick', [], 'ytick', [])
    % xlabel('Occupancy')
    hc = colorbar;
    hc.Label.String = 'Occupancy (sec)';
    
    subplot(m,n,2)
    imagesc(maps2.occ)
    title(names{2})
    set(gca, 'xtick', [], 'ytick', [])
    % xlabel('Occupancy')
    hc = colorbar;
    hc.Label.String = 'Occupancy (sec)';
    
    
    % plot the Tuning curves side by side with the
    for ii = 1:length(overMap.label)
        
        if ip >= ((n*m)/2)-2
            figure(202+ii)
            ip = 1;
        else
            ip = ip+1;
        end
        
        disp(ip)
        subplot(m, n, s1_idx(ip))
        imagesc(maps1.tc{ii})
        axis off
        title(overMap.label{ii})
        hc = colorbar;
        hc.Label.String = 'Rate (spike/sec)';
        
        subplot(m, n, s2_idx(ip))
        imagesc(maps2.tc{ii})
        axis off
        title(["S corr2: " num2str(overMap.corr(ii))])
        hc = colorbar;
        hc.Label.String = 'Rate (spike/sec)';
        
        
    end
    
end




