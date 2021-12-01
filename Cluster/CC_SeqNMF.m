function CC_SeqNMF(ms, behav, k, l)



if nargin < 3
    k = 2; 
    l = 10; 
elseif nargin < 4
    l = 10;
end


%% plug it into SeqNMF
% addpath(PARAMS.code_seqnmf_dir)

Fs =mode(diff(ms.time));

% data_in = ms.FiltTraces';

if ~isfield('ms', 'Binary')
    ms = MS_msExtractBinary_detrendTraces(ms, 2);
end
data_in = ms.Binary';

behav = MS_align_data(behav, ms);


pos(:,1) = interp1(behav.time,behav.position(:,1),ms.time);
pos(:,2) = interp1(behav.time,behav.position(:,2),ms.time);
velocity = interp1(behav.time, behav.speed, ms.time); 

%% limit to laps
movement_thresh = 2.5; % in cm/s
movement_idx = velocity >movement_thresh; % get times when the animal was moving.

left_idx = MS_get_direction(pos(:,1), -.1); % use -threshold for leftbound and + for right.
right_idx = MS_get_direction(pos(:,1), .1); % use -threshold for leftbound and + for right.

left_idx = left_idx & movement_idx; % only keep the indices when they are moving to the left.
right_idx = right_idx & movement_idx;

[L_laps, L_lap_start_idx, L_lap_end_idx] = MS_get_laps(left_idx, floor(1.5*(Fs)),floor(10*(Fs)));
[R_laps, R_lap_start_idx, R_lap_end_idx] = MS_get_laps(right_idx, floor(1.5*(Fs)),floor(10*(Fs)));

% make L and R lap only data sets. 
L_data = data_in(:,L_laps > 0);
R_data = data_in(:,R_laps > 0);

%% for plotting get the position during laps.
L_laps_pos = pos(L_laps > 0,:);
L_laps_time = ms.time(L_laps >0,:);

R_laps_pos = pos(R_laps > 0,:);
R_laps_time = ms.time(R_laps >0,:);

figure(111)
hold on
% plot(ms.time, pos, 'color', [.2 .2 .2])
plot(L_laps_time, L_laps_pos, '.b')
plot(R_laps_time, R_laps_pos, '.r')


%% Fit with seqNMF: most of this is straight out of Mackevicius et al. 2019 https://elifesciences.org/articles/38471
% this_data = data_in;  % which data to use. 
% this_data = R_data;
this_data = L_data;
this_pos = L_laps_pos;

PosFs = Fs; % b/c we already interpolated the position data to the ms.time. 

% break data into training set and test set
splitN = floor(size(this_data,2)*.75); 
splitS = floor(size(this_pos,2)*.75);

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
K = 2;
L = 10;
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


end