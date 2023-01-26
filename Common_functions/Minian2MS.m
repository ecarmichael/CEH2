function ms = Minian2MS(data_dir, plot_flag, save_flag)
%% Minian2Mat: Read the netCDF output from Minian. based on function found here: https://groups.google.com/g/miniscope/c/fW7xGqWLd4E
%
%
%
%    Inputs:
%    - data_dir: [string]  path to data
%
%
%
%    Outputs:
%    - ms: [struct] contains the outputs from Minian
%
%
%
%
% EC 2022-12-18   initial version
%
%
%
%% initialize

if nargin <1
    data_dir = cd;
    plot_flag = 1;
    save_flag = 0; 
elseif nargin <2
    plot_flag = 1;
    save_flag = 0;
elseif nargin < 3
    save_flag = 0;
end

if isempty(dir('*.nc'))
    error('No .nc file detected in current directory')
end

%% load everything
fname = [data_dir filesep 'minian_dataset.nc'];
ms = [];

ms.SFP  = ncread(fname,'A');
%compute spf centroids
for ii = size(ms.SFP, 3):-1:1
    ilabel = logical(imfill((ms.SFP(:,:,ii)), 'holes'));
    c_out = regionprops(ilabel,'centroid');
    
    ms.centroids(:,ii) = c_out.Centroid;
end

ms.Deconv  = ncread(fname,'C');
ms.Spikes  = ncread(fname,'S');
ms.c0 = ncread(fname,'c0'); % convolved traces?
ms.background = ncread(fname,'b');
ms.b0 = ncread(fname,'b0');
ms.background_trace = ncread(fname,'f');
ms.RawTraces = ncread(fname,'YrA');

ms.frame = double(ncread(fname,'frame'));

ms.units = double(ncread(fname,'unit_id'));
ms.unit_labels = double(ncread(fname,'unit_labels'));
ms.unit_labels = ms.unit_labels+1;

ms.max_proj = ncread(fname,'max_proj');


%Less usefull variables: height, width, frame, motion
ms.height = double(ncread(fname,'height'));
ms.width = double(ncread(fname,'width'));
ms.motion = double(ncread(fname,'motion'));

% load the timestamps. 
TS = readtable('timeStamps.csv');

tvec = table2array(TS(:,2));

% correct for offsets if needed
if tvec(1) ~= 0
    tvec = tvec+abs(tvec(1));
end
ms.tvec = tvec./1000; % convert to seconds

%Variables; animal, session, shift_dim are of an unsupported datatype and
% cannot be extracted with ncread


% %% Make version of variables with bad cells (label =-1) removed
% indx = 0 <= ms.unit_labels;
%
% ids = ms.unit_labels(indx);
% A = ms.SFP(:,:,indx);
% C = ms.RawTraces(:,indx);
% S = ms.Sdata(:,indx);
% c0 = ms.c0data(:,indx);
% YrA = ms.YrAdata(:,indx);

%% quick plot to check data

if plot_flag
    figure(1919)
    clf
    if length(ms.units) < 20
        nCells  = length(ms.units);
    else
        nCells = 20;
    end
    c_ord = linspecer(nCells);
    
    subplot(2,2,2)
    imagesc(mean(ms.SFP(:,:,:),3))
    hold on

    for ii = 1:nCells
        [row, col] = find(ms.SFP(:,:,ii) == max(max(ms.SFP(:,:,ii))));
        scatter(col, row,50, 'MarkerEdgeColor', c_ord(ii,:), 'LineWidth', 2)
    end
    title(['nCells: ' num2str(length(ms.units))]);
    xlim([min(ms.centroids(1,:))-min(ms.centroids(1,:))*.2  max(ms.centroids(1,:))+max(ms.centroids(1,:))*.2])
    ylim([min(ms.centroids(2,:))-min(ms.centroids(2,:))*.2  max(ms.centroids(2,:))+max(ms.centroids(2,:))*.2])

    subplot(2,2,4)
    cla
    all_cord = parula(floor(length(ms.units)*1.2));%zeros(length(ms.units),3);%repmat(linspecer(5), 100, 1); 
    hold on
    mult_fac = 5; 
    for ii = 1:length(ms.units)
        plot(ms.tvec, zscore(ms.Deconv(:,ii))+ii*mult_fac,  'color', all_cord(ii,:), 'linewidth', 1.5)

%         plot(ms.frame*(1/30), zscore(ms.Deconv(:,ii))+ii*10,'color', all_cord(ii,:), 'linewidth', 1)
%         plot(ms.frame*(1/30), (ms.Spikes(:,ii)*5)+ii*10,'color', 'k')
        
    end
    title('zscored deconvolved traces (all cells)')
    ylabel('cell ID')
    xlabel('time (s)')
    set(gca,'ytick', 0:100:length(ms.units)*mult_fac, 'YTickLabel', (0:100:length(ms.units)*mult_fac)/mult_fac, 'TickDir', 'out')
    ylim([mult_fac (length(ms.units)+2)*mult_fac])
    xlim([ms.tvec(1) ms.tvec(end)])
    
    
      subplot(2,2,[1 3])
    hold on
    for ii = 1:nCells
        plot(ms.tvec, ms.RawTraces(:,ii)+ii*10,'color', [0.8 0.8 0.8], 'linewidth', 1)
        plot(ms.tvec, ms.Deconv(:,ii)+ii*10,'color', c_ord(ii,:), 'linewidth', 2)
        plot(ms.tvec, (ms.Spikes(:,ii)*5)+ii*10,'color', 'k')
        
    end
    title('Denoise and deconvolved traces for sample cells')
    ylabel('cell ID')
    xlabel('time (s)')
    set(gca,'ytick', 0:10:nCells*10, 'YTickLabel', 0:nCells+1, 'TickDir', 'out')
    ylim([10 (nCells+2)*10])
    xlim([ms.tvec(1) ms.tvec(end)])
    
%     subplot(2,2,4)
%     hold on
%     for ii = 1:nCells
%         plot(ms.frame, ms.c0(:,ii)+ii,'color', c_ord(ii,:))
%     end
%     title('Convolved traces for sample cells')
%     xlabel('cell ID')
%     ylabel('time (s)')
maximize
saveas(gcf, 'Minian_Screener.png');
saveas(gcf, 'Minian_Screener.fig')

end

if save_flag
    save('minian_ms', 'ms')
end










end