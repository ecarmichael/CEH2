function [A_temp, A_proj, data_h, tvec] =MS_PCA_ICA_only(ms, move_idx, binsize, method, opts)

if nargin < 3
    binsize = .5;
    method = [];
    opts = [];
elseif nargin < 4
    method = [];    
    opts = []; 
elseif nargin < 5
    opts  = []; 
end

if isempty(method)
    method = 'binary';
end

if isempty(opts)
    opts.threshold.method = 'MarcenkoPastur';
    opts.Patterns.method = 'ICA';
    opts.Patterns.number_of_iterations = 500;
end
    

rng(123, 'twister')

%%

%% follow grosmark et al. method of deconv preprocessing
if strcmp(method, 'grossmark')
    Csp = ms.deconv./ms.denoise;
    Csp = Csp > 0.01;
    ms.Csp = Csp;
    data = Csp;
else
    data = ms.Binary;
end





% %% bin and convolve
% binsize = 0.1; % in seconds, so everything else should be seconds too
% gauss_window = 1./binsize; % 1 second window
% gauss_SD = 0.5./binsize; % 0.02 seconds (20ms) SD
% gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
% gau_sdf = conv2(Csp,gk,'same'); % convolve with gaussian window
%
% gau_z = zscore(gau_sdf, [], 2);

%% try again with histc
% binsize = .5;
tbin_edges = ms.time(1)/1000:binsize:ms.time(end)/1000; % vector of time bin edges (for histogram)
tbin_centers = tbin_edges(1:end-1)+binsize/2; % vector of time bin centers (for plotting)

data_h = [];
for ii = size(ms.Binary,2):-1:1
    
    %     this_cell = ms.time(find(Csp(: ,ii) & move_idx))/1000;
    this_cell = ms.time(find(ms.Binary(: ,ii) & move_idx))/1000;
    
    spk_count = histc(this_cell,tbin_edges); % get spike counts for each bin
    spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.
    
    data_h(:,ii) = spk_count;
%     data_h(:,ii) = zscore(spk_count);
end

tvec = tbin_centers;


%% try the assembly code....
rng(123, 'twister')

% get the assemblies

A_temp = assembly_patterns(data_h', opts);

if size(A_temp,1) == size(data_h,2)
    disp('Templates match # cells')
else 
    error('Input data has more cells than timepoints');
end


% get the projections
rng(123, 'twister')
A_proj = assembly_activity(A_temp,data_h');


%% shuffle distribution for assemblies
% rng(123,'twister')
% nShuff = 100;
% wake_shuff_mat = [];
% 
% Ass_shuff = NaN(1,nShuff);
% for iS = nShuff:-1:1
%     tic
%     shuff_data = NaN(size(data_h));
%     for ic = 1:size(data_h,2)
%         shuff_data(:,ic) = circshift(data_h(:,ic), floor(MS_randn_range(1,1,1,size(data_h,1))));
%     end
%     
%     this_ass = assembly_patterns(shuff_data');
%     if ~isempty(this_ass)
%         S_prog = assembly_activity(this_ass,shuff_data');
%         
%         wake_shuff_mat(iS,:) =  S_prog(1,:);
%         keep_idx(iS) = 1;
%     else
%         wake_shuff_mat(iS,:) = NaN;
%         keep_idx(iS) = 0;
%     end
%     %     for ii = size(this_ass,2):-1:1
%     
%     if sum(max(this_ass) > 0.2) >0
%         Ass_shuff(iS) = sum(max(this_ass) > 0.2);
%     else
%         Ass_shuff(iS) = 0;
%     end
%     %     end
%     fprintf('Shuff # %.0f found %.0f assemblies and took %2.2f seconds\n', iS, size(this_ass,2), toc)
% end
% 
% W_threshold = prctile(wake_shuff_mat(wake_shuff_mat >0), 99, 'all');

