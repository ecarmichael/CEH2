function ms = MS_minian2ms(minian_fname, plot_flag)
%% MS_minian2ms: extract data from minian output in the netCDF format and exports it in the ms struct akin to CNMFe
%
%
%
%    Inputs: 
%    - minian_fname: [string] name of the file to load
%
%    - plot_flag: logical 
%
%    Outputs: 
%    - ms: [struct]  data structure with subfields for time, traces,
%    spikes, and footprints. 
%
%
%
%
% EC 2022-04-18   initial version 
%
%
%
%% initialize
if nargin <2
    plot_flag =0; 
end
    

fprintf('<strong>%s</strong>: loading data file: <strong>%s</strong>...\n', mfilename, minian_fname); 


A = ncread(minian_fname,'A');                 
C = ncread(minian_fname,'C');
S = ncread(minian_fname,'S');

ms.dirname = cd; % where did you find the file?
ms.numFrames = size(C); % get the number of frames in the data. 
ms.maxFramesPerFile = 1000; % default for the miniscope software. 
ms.height = size(A, 2); 
ms.width = size(A,1); 
ms.shifts = ncread(minian_fname,'motion')'; % get the shifts
ms.meanFrame = mean(A,3)'; % compute the mean frame. 

fprintf('<strong>%s</strong>: computing centroids...\n', mfilename); 
for ii = size(A,3):-1:1
    cent = regionprops(true(size(A(:,:,ii))), A(:,:,ii),  'WeightedCentroid'); 
        ms.Centroids(ii,:) = cent.WeightedCentroid; 
end

ms.SPFs = A; 
% get frames and correct for python indexing. 
ms.frameNum = ncread(minian_fname,'frame');
ms.frameNum = ms.frameNum+1; %compensate for python indexing
max_proj = ncread(minian_fname, 'max_proj'); 



%% plot a check if needed. 
if plot_flag
    c_ord = winter(size(C, 2)); 
    s_ord = autumn(size(C, 2)); 
    
    figure(99999)
    clf
    
    subplot(1,5,1:2)
%     imagesc(1:ms.w ms.max_proj'); hold on
    set(gca, 'YDir', 'normal'); 
    
    subplot(1,5,3:5)
    hold on
    for iC = 1:5:size(C,2)
       plot(ms.frameNum, (C(:,iC)*.1)+iC, 'color', c_ord(iC,:))
       nan_idx = S(:,iC) <=0; 
       plot(ms.frameNum(~nan_idx), (S(~nan_idx,iC)*.1)+iC,'.k')

        
    end

    
    
end
