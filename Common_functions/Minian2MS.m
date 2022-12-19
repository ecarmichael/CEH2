function ms = Minian2MS(data_dir)
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
end

if ~isempty(dir('*.nc'))
    error('No .nc file detected in current directory')
end

%% load everything
fname = [data_dir filesep 'minian_dataset.nc']; 
ms = [];

ms.SFP  = ncread(fname,'A');                  
ms.Deconv  = ncread(fname,'C');
ms.Spikes  = ncread(fname,'S');
ms.c0 = ncread(fname,'c0');
ms.background = ncread(fname,'b');
ms.b0 = ncread(fname,'b0');
ms.background_trace = ncread(fname,'f');
% RawTraces = ncread(fname,'YrA');

ms.frame = ncread(fname,'frame');

ms.units = ncread(fname,'unit_id');
ms.unit_labels = ncread(fname,'unit_labels');
ms.unit_labels = ms.unit_labels+1; 

ms.max_proj = ncread(fname,'max_proj');


%Less usefull variables: height, width, frame, motion
ms.height = ncread(fname,'height');
ms.width = ncread(fname,'width');
ms.motion = ncread(fname,'motion');

%Variables; animal, session, shift_dim are of an unsupported datatype and
% cannot be extracted with ncread


%% Make version of variables with bad cells (label =-1) removed
indx = 0 <= ms.unit_labels;

ids = unit_labels(indx);
A = ms.SFP(:,:,indx);
C = ms.RawTraces(:,indx);
S = ms.Sdata(:,indx);
c0 = ms.c0data(:,indx);
% YrA = ms.YrAdata(:,indx);

%% quick plot to check data

if plot_flag 
    figure(1919)
    clf
        c_ord = linspecer(10); 

    subplot(2,2,2)
    imagesc(mean(ms.SFP(:,:,:),3))
    hold on
    for ii = 1:10
        [row, col] = find(ms.SFP(:,:,ii) == max(max(ms.SFP(:,:,ii)))); 
        scatter(col, row,50, 'MarkerEdgeColor', c_ord(ii,:), 'LineWidth', 2)
    end
    
    subplot(2,2,[1 3])
    hold on
    for ii = 1:10
        plot(ms.frame, ms.Deconv(:,ii)+ii*10,'color', c_ord(ii,:), 'linewidth', 2)
        plot(ms.frame, (ms.Spikes(:,ii)*5)+ii*10,'color', 'k')

    end
    set(gca,'ytick', 0:10:100, 'YTickLabel', 0:11)
    
     subplot(2,2,4)
    hold on
    c_ord = linspecer(10); 
    for ii = 1:10
        plot(ms.frame, ms.c0(:,ii)+ii*10,'color', c_ord(ii,:))
    end
end
    
    
    
    
    
    
    
    
    
    
end