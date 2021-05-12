function h = MS_plot_spatial_cell(SI_in, cell_id)
%% MS_plot_spatial_cell: plot the spatial properities of a cell
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


M = 3; % rows
N = 4; % columns

 figure

c_num = 0; % track the cell number for plotting across figures.

% loop over cells.
for iC = cell_id
    c_num = c_num +1; % for plotting. 
    
    % make a new figure if beyong 12. 
    if c_num <= N*M
       c_num = 1;
       figure; % make a new figure; 
        
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
    subplot(M, N, c_num) % N*4+4:N*4+6
    imagesc(this_cell.spatial.place.Sig_map);
    set(gca, 'YDir', 'normal');
    xlabel('position (cm)');
    ylabel('position (cm)');
    tmp=get(gca,'position');
    set(gca,'position',[tmp(1) tmp(2) tmp(3) (X_bin_centers(end) /Y_bin_centers(end))*tmp(4)])
    
    
    % plot the MI and p value for the cell.
    subplot(M, N, N+5)
    text(0, 1*max(ylim), 'Place', 'HorizontalAlignment', 'left', 'color', 'K', 'fontweight', 'bold')
    text(0, .8*max(ylim), {'MI:'; num2str(SI.place.MI(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
    text(0, .4*max(ylim), {'split corr:'; num2str(SI.place.split.Stability_corr(iC),3)}, 'HorizontalAlignment', 'left', 'color', 'K')
    axis off
    
end % cells.

