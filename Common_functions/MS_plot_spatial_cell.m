function h = MS_plot_spatial_cell(SI_in, cell_id)
%% MS_plot_spatial_cell: plot the spatial properities across cells
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
% EC 2021-05-05   initial version
%
%% initialize
if nargin <2
    cell_id = 1:length(SI_in);
end
%%  basic fig properties

% subplot properties
M = 3; % rows
N = 3; % columns

figure
set(gcf, 'position', [680 111 1080 864]); % a nice square-ish plot
c_num = 0; % track the cell number for plotting across figures.

% get the spatial bins
X_bin_centers = SI_in(1).cfg.X_bins +  SI_in(1).cfg.X_bins/2;
X_bin_centers = X_bin_centers(1:end-1); 

Y_bin_centers = SI_in(1).cfg.Y_bins +  SI_in(1).cfg.Y_bins/2;
Y_bin_centers = Y_bin_centers(1:end-1); 

% loop over cells.
for iC = cell_id
    c_num = c_num +1; % for plotting. 
    
    % make a new figure if beyong 12. 
    if c_num >= N*M
       c_num = 1;
       figure; % make a new figure; 
        set(gcf, 'position', [680 111 1080 864])
    end
    % run this cell
    this_cell = SI_in(iC);
 
    % add in the 2D place/spatial information?
    subplot(M, N, c_num) % N*4+4:N*4+6
    if max(X_bin_centers) < max(Y_bin_centers)
        imagesc(this_cell.cfg.X_bins, this_cell.cfg.Y_bins, this_cell.spatial.place.Sig_map);
    else
            imagesc(this_cell.cfg.Y_bins, this_cell.cfg.X_bins, this_cell.spatial.place.Sig_map');
    end
    set(gca, 'YDir', 'normal');
    xlabel('position (cm)');
    ylabel('position (cm)');
%     tmp=get(gca,'position');
%     set(gca,'position',[tmp(1) tmp(2) tmp(3) (X_bin_centers(end) /Y_bin_centers(end))*tmp(4)])
    
    title(['Cell: ' num2str(iC) ' MI: ' num2str(this_cell.spatial.place.MI,2) ' split cor: ' num2str(this_cell.spatial.place.split.Stability_corr,2)], 'fontsize', 10)
    
    % add some session info
    if c_num == 1
        text(-50*abs(min(xlim)),1.2*max(ylim),['Subject: ' this_cell.finfo.subject ' | Date: ' this_cell.finfo.date ' '  this_cell.finfo.task],...
            'fontweight', 'bold', 'fontname', 'helvetica', 'fontsize', 12);
    end
    
end % cells.

