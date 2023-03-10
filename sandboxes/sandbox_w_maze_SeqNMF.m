%% sandbox_w_maze_SeqNMF
addpath('C:\Users\ecarm\Documents\GitHub\seqNMF')

%% load some data

load('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1252\2021_12_16_pv1252_MZD3\ms_trk.mat');

%% DLC conversion

pos = MS_DLC2TSD('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\Raw\2021_12_16_pv1252_MZD3\H12_M50_S24_MAZE');

pos = MS_align_data_uni(pos, ms_trk);
figure(1)
clf
hold on
plot(pos.data(1,:), pos.data(2,:), '.')

%%

ms = [];
ms.time = ms_trk.time;
ms.Binary  = ms_trk.Binary; 
ms.deconv  = ms_trk.deconv; 
ms.denoise  = ms_trk.denoise; 
clear ms_trk

Csp = ms.deconv./ms.denoise; 
Csp = Csp > 0.01; 

% bin and convolve
binsize = 0.1; % in seconds, so everything else should be seconds too
gauss_window = 1./binsize; % 1 second window
gauss_SD = 0.5./binsize; % 0.02 seconds (20ms) SD
gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
gau_sdf = conv2(Csp,gk,'same'); % convolve with gaussian window

gau_z = zscore(gau_sdf, [], 2); 

gau_z(gau_z <0) = 0; 

ms.gau_z = gau_z; 
clear gau_z Csp
%% plot the gaussian smoothed Csp
figure(101)
ax(1)= subplot(5,1,1);
hold on
plot(ms.time, pos.data(1,:)); 
plot(ms.time, pos.data(2,:)); 
xlim([ms.time(1) ms.time(end)])

ax(2) = subplot(5,1,2:5);
imagesc(ms.time, 1:size(Csp,2), gau_z')

linkaxes(ax, 'x')

    
    
clear ms_trk; 



%% Fit with seqNMF: most of this is straight out of Mackevicius et al. 2019 https://elifesciences.org/articles/38471

Fs = mode(diff(ms.time));

this_data = ms.gau_z'; 
this_pos = pos_i.data(1,:);

PosFs = Fs; % b/c we already interpolated the position data to the ms.time. 

% break data into training set and test set
splitN = floor(size(this_data,2)*.70); 
splitS = floor(size(this_pos,2)*.70);

% make the test/trian sets.  
trainNEURAL = this_data(:,1:splitN); 
trainPOS = this_pos(:,1:splitS); 
testNEURAL = this_data(:,(splitN+1):end); 
testPOS = this_pos(:,(splitS+1):end); 
%% plot one example factorization
rng(235); % fixed rng seed for reproduceability
X = trainNEURAL;
figure
% loop L values
% Ls = fliplr([0.5 1 2 10 50 100 1000]);
% 
% for iS = length(Ls):-1:1
K = 3;
L = 15;
% L = Ls(iS); % units of seconds
Lneural = ceil(L*Fs);  
% Lsong = ceil(L*SONGfs);
shg
subplot(2,2,1)
display('Running seqNMF...')
[W, H, ~,loadings,power]= seqNMF(X,'K',K,'L',Lneural,...
            'lambdaL1W', .1, 'lambda', .005, 'maxiter', 100, 'showPlot', 1,...
            'lambdaOrthoW', 0); 
        
        %%
p = .05; % desired p value for factors

display('Testing significance of factors on held-out data')
[pvals,is_significant] = test_significance(testNEURAL,W,p);

W = W(:,is_significant,:); 
H = H(is_significant,:); 
fprintf('Found %d/%d significant factors\n', sum(is_significant), length(is_significant))
% all_sweeps{iS}.K = K;
% all_sweeps{iS}.L = L;
% all_sweeps{iS}.W = W;
% all_sweeps{iS}.H = H;
% all_sweeps{iS}.sig = is_significant;
% all_sweeps{iS}.Train = trainNEURAL;
% all_sweeps{iS}.Trest = testNEURAL;
% 
% saveas(gcf,['Seq_Sweeps' filesep 'Sweep_' num2str(L)], 'png')
close all
% clearvars -except all_sweeps trainNEURAL testNEURAL Ls

% end
%%
% plot, sorting neurons by latency within each factor
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
indSort = hybrid(:,3);
tstart = 1; % plot data starting at this timebin
figure;
WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), ...
    0,trainPOS(:,floor(tstart*Fs/PosFs):end))
title('Significant seqNMF factors, with raw data')

figure;
WHPlot(W(indSort,:,:),H(:,tstart:end), ...
    helper.reconstruct(W(indSort,:,:),H(:,tstart:end)),...
    0,trainPOS(:,floor(tstart*Fs/PosFs):end))
title('SeqNMF reconstruction')