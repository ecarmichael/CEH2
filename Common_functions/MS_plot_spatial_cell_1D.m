function h = MS_plot_spatial_cell_1D(SI_in, cell_id)
%% MS_plot_spatial_cell_1D: plot the spatial properities of a cell in 1D
%
%
%
%    Inputs:
%    -SI_in: [struct]  output from Spatial_screener_info.m
%
%    - cell_id: [1 x N vector] the cell index to process. can be single or
%    a vector.
%
%    Outputs:
%    - h [handle]
%
%
%
%
% EC 2021-05-05   initial version
%
%
%
%% initialize
if nargin <2
    cell_id = 1:length(SI_in);
end
%%  basic fig properties

M = 12; % rows
N = 4; % columns
sub_idx = reshape(1:M*N, N,M)'; 


figure
set(gcf, 'position', [232   110   1000   840]);

c_num = 0; % track the cell number for plotting across figures.
% loop over cells.
for iC = cell_id
    c_num = c_num +2; % for plotting. 
    
    % make a new figure if beyong 12. 
    if c_num <= floor(M/2)
       c_num = 2;
       figure; % make a new figure; 
       set(gcf, 'position', [232   110   1000   840]); 
    end
    % run this cell
    this_cell = SI_in(iC);
    
    %%% title information
%     subplot(M, N, 4:5) % title information. Top right corner,
%     ylim([0 10])
%     text(0,8,['Cell id: ' num2str(iC)], 'fontweight', 'bold')
%     text(0,6,['Subject: ' this_cell.fname.subject]);
%     text(0,4,['Session: ' this_cell.fname.task]);
%     text(0,2,['Date: ' this_cell.fname.date]);
    
%     axis off
 
    % add in the 2D place/spatial information?
    this_sub = sub_idx(c_num-1:c_num,1:3); 
    subplot(M, N, this_sub(:)'); 
    imagesc(this_cell.cfg.X_bins,1,  this_cell.spatial.place.posterior');
    set(gca, 'YDir', 'normal', 'ytick', []);
    xlabel('position (cm)');
    colormap([0,0,0; parula(128)]); 
    text(0, 1.2*max(ylim), ['Place | MI: ' num2str(this_cell.spatial.place.MI,3) ' | split corr: ' num2str(this_cell.spatial.place.split.Stability_corr,3)], 'HorizontalAlignment', 'left', 'fontweight', 'bold', 'color', 'K')
    colorbar;
    
    % split half 1
        this_sub = sub_idx(c_num-1,4); 
    subplot(M, N, this_sub(:)'); 
    imagesc(this_cell.cfg.X_bins,1,  this_cell.spatial.place.split.S1_Sig_map');
    set(gca, 'YDir', 'normal', 'ytick', []);
    xlabel('position (cm)');
    colormap([0,0,0; parula(128)]); 
    text(0, 1.2*max(ylim), ['MI: ' num2str(this_cell.spatial.place.split.S1_MI,3)], 'HorizontalAlignment', 'left', 'fontweight', 'bold', 'color', 'K')

        % split half 2
        this_sub = sub_idx(c_num,4); 
    subplot(M, N, this_sub(:)'); 
    imagesc(this_cell.cfg.X_bins,1,  this_cell.spatial.place.split.S2_Sig_map');
    set(gca, 'YDir', 'normal', 'ytick', []);
    xlabel('position (cm)');
    colormap([0,0,0; parula(128)]); 
    text(0, 1.2*max(ylim), ['MI: ' num2str(this_cell.spatial.place.split.S2_MI,3)], 'HorizontalAlignment', 'left', 'fontweight', 'bold', 'color', 'K')

    
    
end % cells.

