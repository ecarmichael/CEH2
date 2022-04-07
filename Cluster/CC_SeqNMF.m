function CC_SeqNMF(ms, behav, k, l, lambda, save_dir )



if nargin < 3
    k = 2;
    l = 10;
    lambda = 0.005;
    save_dir = [cd filesep 'Seq_out_L' strrep(num2str(l), '.', 'p')];
elseif nargin < 4
    l = 10;
    lambda = 0.005;
    save_dir = [cd filesep 'Seq_out_L' strrep(num2str(l), '.', 'p')];
elseif nargin < 6
    lambda = 0.005;
    save_dir = [cd filesep 'Seq_out_L' strrep(num2str(l), '.', 'p')];
    elseif nargin < 7
    save_dir = [cd filesep 'Seq_out_L' strrep(num2str(l), '.', 'p')];
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

[L_laps, ~, ~] = MS_get_laps(left_idx, floor(1.5*(Fs)),floor(10*(Fs)));
[R_laps, ~, ~] = MS_get_laps(right_idx, floor(1.5*(Fs)),floor(10*(Fs)));

% make L and R lap only data sets. 
L_data = data_in(:,L_laps > 0);
R_data = data_in(:,R_laps > 0);

%% for plotting get the position during laps.
L_laps_pos = pos(L_laps > 0,:);
L_laps_time = ms.time(L_laps >0,:);
L_laps_velo = velocity(L_laps >0,:);

R_laps_pos = pos(R_laps > 0,:);
R_laps_time = ms.time(R_laps >0,:);
R_laps_velo = velocity(R_laps >0,:);

figure(111)
clf
subplot(2,1,1)
hold on
% plot(ms.time, pos, 'color', [.2 .2 .2])
plot(L_laps_time, L_laps_pos(:,1), '.b')
plot(R_laps_time, R_laps_pos(:,1), '.r')
legend({'left'; 'right'}, 'Location', 'northeast', 'Orientation', 'horizontal')
subplot(2,1,2)
% h0 = axes;
hold on
scatter(ms.time, behav.position(:,1), 5,'k');
scatter(ms.time(movement_idx), behav.position(movement_idx,1), 5,velocity(movement_idx));
colorbar('Location', 'northoutside')
colormap('parula')
ylim([min(behav.position(:,1)) max(behav.position(:,1))*1.3])
% set(gca, 'xtick', [], 'ytick', [])
% h1 = axes;
% scatter(ms.time(left_idx), behav.position(left_idx,1), 10, velocity(left_idx), 'filled');
% h2 = axes; 
% scatter(ms.time(right_idx), behav.position(right_idx,1), 10, velocity(right_idx), 'filled');
% linkaxes([h0, h1,h2])
% h2.Visible = 'off'; h2.XTick = []; h2.YTick = []; h2.XTickLabel = []; h2.YTickLabel = []; 
% h1.Visible = 'off'; h1.XTick = []; h1.YTick = []; h1.XTickLabel = []; h1.YTickLabel = []; 
% h0.Visible = 'off'; h0.XTick = []; h0.YTick = []; h0.XTickLabel = []; h0.YTickLabel = []; 

% colormap(h1, 'winter')
% colormap(h2, 'hot')
% hold on
% scatter(ms.time(~left_idx | ~right_idx), behav.position(~left_idx | ~right_idx,1), 5, velocity(~left_idx | ~right_idx), 'k');


%% Fit with seqNMF: most of this is straight out of Mackevicius et al. 2019 https://elifesciences.org/articles/38471
% this_data = data_in;  % which data to use. 
% this_data = R_data;
% fprintf('Running Left laps\n')
% this_data = L_data;
% this_pos = L_laps_pos;
% this_velo = L_laps_velo; 

fprintf('Running all laps\n')
% this_data = R_data;
% this_pos = R_laps_pos;
% this_velo = R_laps_velo; 

fprintf('Running all laps\n')
this_data = data_in(:,movement_idx);
this_pos = ms.time(movement_idx);
this_velo = velocity(movement_idx); 

PosFs = Fs; % b/c we already interpolated the position data to the ms.time. 

% break data into training set and test set
splitN = floor(size(this_data,2)*.75); 
splitS = floor(size(this_pos,1)*.75);

% convert 1D position into maxtix
bin_s = 2.5; pos_bins = 0:bin_s:100; 
pos_mat = zeros(length(this_pos),length(pos_bins));
for ii = length(this_pos):-1:1
    pos_mat(ii,nearest_idx3(this_pos(ii,1), pos_bins)) = this_velo(ii); 
end
pos_mat = (pos_mat); 
% make the test/trian sets.  
trainNEURAL = this_data(:,1:splitN); 
trainPOS = pos_mat(1:splitS,:); 
testNEURAL = this_data(:,(splitN+1):end); 
testPOS = pos_mat((splitS+1):end,:); 


%% plot one example factorization
rng(235); % fixed rng seed for reproduceability
X = trainNEURAL;

% loop L values
% Ls = fliplr([0.5 1 2 10 50 100 1000]);
% 
% for iS = length(Ls):-1:1
% K = 2;
% L = 10;
% L = Ls(iS); % units of seconds
Lneural = ceil(l*Fs);  
% Lsong = ceil(L*SONGfs);


%% find lambda
if isempty(lambda)
nLambdas = 20; % increase if you're patient
K = k; 
X = trainNEURAL;
lambdas = sort([logspace(-1,-5,nLambdas)], 'ascend'); 
loadings = [];
regularization = []; 
cost = []; 
for li = 1:length(lambdas)
    [N,T] = size(X);
    [W, H, ~,loadings(li,:),power]= seqNMF(X,'K',K,'L',Lneural,...
        'lambdaL1W', .1, 'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0); 
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X,W,H);
    display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
end

windowSize = 3; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
Rs = filtfilt(b,a,regularization); 
minRs = prctile(regularization,10); maxRs= prctile(regularization,90);
Rs = (Rs-minRs)/(maxRs-minRs); 
R = (regularization-minRs)/(maxRs-minRs); 
Cs = filtfilt(b,a,cost); 
minCs =  prctile(cost,10); maxCs =  prctile(cost,90); 
Cs = (Cs -minCs)/(maxCs-minCs); 
C = (cost -minCs)/(maxCs-minCs); 

clf; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
[~, opt_l_idx] = min(abs(C - R)); 
opt_lambda = lambdas(opt_l_idx);
plot(lambdas(opt_l_idx), C(opt_l_idx), 'xk', 'MarkerSize', 10)
else
   opt_lambda =  lambda; 
    
end
%% run Seq
% shg
% subplot(2,2,1)
% fprintf('Running seqNMF (k = %.0f ; l = %.0fsec)...\n', k, l)
% [W, H, ~,loadings,power]= seqNMF(X,'K',k,'L',Lneural,...
%             'lambdaL1W', .1, 'lambda', opt_lambda, 'maxiter', 100, 'showPlot', 1,...
%             'lambdaOrthoW', 0); 
        
        %%
% p = .05; % desired p value for factors
% 
% display('Testing significance of factors on held-out data')
% [pvals,is_significant] = test_significance(testNEURAL,W,p);
% 
% W = W(:,is_significant,:); 
% H = H(is_significant,:); 
% fprintf('Found %d/%d significant factors\n', sum(is_significant), length(is_significant))
% % all_sweeps{iS}.K = K;
% % all_sweeps{iS}.L = L;
% % all_sweeps{iS}.W = W;
% % all_sweeps{iS}.H = H;
% % all_sweeps{iS}.sig = is_significant;
% % all_sweeps{iS}.Train = trainNEURAL;
% % all_sweeps{iS}.Trest = testNEURAL;
% % 
% % saveas(gcf,['Seq_Sweeps' filesep 'Sweep_' num2str(L)], 'png')
% close all
% % clearvars -except all_sweeps trainNEURAL testNEURAL Ls

% end
%% sort and reconstruc 
% % plot, sorting neurons by latency within each factor
% [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
% indSort = hybrid(:,3);
% tstart = 1; % plot data starting at this timebin
% figure;
% WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), ...
%     1,trainPOS(:,floor(tstart*Fs/PosFs):end)')
% title('Significant seqNMF factors, with raw data')
% 
% figure;
% WHPlot(W(indSort,:,:),H(:,tstart:end), ...
%     helper.reconstruct(W(indSort,:,:),H(:,tstart:end)),...
%     1,trainPOS(:,floor(tstart*Fs/PosFs):end)')
% title('SeqNMF reconstruction')
% 


%% testing

% figure(111)
% clf
% 
% ax(1) = subplot(5,1,1);
% plot((H(:,tstart:end))')
% xlim([1 length(H(:,tstart:end))'])
% 
% ax(2) = subplot(5,1,2);
% imagesc(X(indSort,tstart:end));
% 
% ax(3) = subplot(5,1,3);
% imagesc((trainPOS(:,floor(tstart*Fs/PosFs):end)'))
% 
% linkaxes(ax, 'x')


%% interative parts based 
nIter = 20; % increase if patient

lambdaOrthoH = 0;  
lambdaOrthoW = 1; % favor parts-based (these can take any value, don't need to be zero and one)


fprintf('Running seqNMF multiple times for lambda=%0.4f\n', opt_lambda)
Ws = []; Hs = []; loadings = []; pvals = []; is_significant = []; hybrid = []; 

for iteri = nIter:-1:1
    tic
 [W, H, ~,loadings(iteri,:),power]= seqNMF(X,'K',k,'L',Lneural,...
            'lambdaL1W', .1, 'lambda', opt_lambda,'lambdaOrthoH', lambdaOrthoH, 'lambdaOrthoW', lambdaOrthoW, 'maxiter', 100, 'showPlot', 0); 
    p = .05;
    [pvals(iteri,:),is_significant(iteri,:)] = test_significance(testNEURAL,W,p);
    Ws{iteri} = W(:,is_significant(iteri,:)==1,:); 
    Hs{iteri} = H(is_significant(iteri,:)==1,:); 
    if sum(is_significant(iteri,:)) > 0
        [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
        indSort = hybrid(:,3);
        tstart = 1;
        figure(iteri); clf; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 1,trainPOS(:,floor(tstart*Fs/PosFs):end)')
    else
        hybrids = [];
    end
    display(['seqNMF run ' num2str(iteri) '/' num2str(nIter)])
    hybrids{iteri} = hybrid; 
    fprintf('Iter # .0f ', iteri);  
    toc
end

figure(1010); hold on
h = histogram(sum(is_significant,2), 'edgecolor', 'w', 'facecolor', .7*[1 1 1]); 
h.BinCounts = h.BinCounts/sum(h.BinCounts)*100; 
xlim([0 10]); 
xlabel('# significant factors')
ylabel('% seqNMF runs')

%% collect outputs
Seq_out.k = k;
Seq_out.l = l;
Seq_out.lambda = opt_lambda; 
Seq_out.lambdaOrthoH = lambdaOrthoH;
Seq_out.lambdaOrthoW = lambdaOrthoW; 
Seq_out.this_data = this_data; 
Seq_out.this_pos = this_pos;
Seq_out.this_velo = this_velo;
Seq_out.trainNEURAL = trainNEURAL; 
Seq_out.trainPOS= trainPOS; 
Seq_out.testNEURAL = testNEURAL; 
Seq_out.testPOS = testPOS; 
Seq_out.W = Ws;
Seq_out.H = Hs;
Seq_out.hybrid = hybrid; 
Seq_out.loadings = loadings;
Seq_out.is_significant = is_significant;
Seq_out.pvals = pvals; 

%% save come ouputs
mkdir(save_dir)
save([save_dir filesep 'Seq_out'], 'Seq_out',  '-v7.3'); 

figure(1010)
saveas(gcf, [save_dir filesep 'Sig_seqs_bar'], 'fig')
pause(1)
saveas(gcf, [save_dir filesep 'Sig_seqs_bar'], 'png')

figure(111)
saveas(gcf, [save_dir filesep 'position'], 'fig')
pause(1)
saveas(gcf, [save_dir filesep 'position'], 'png')

count = 0; 
for ii = 1:10
    if count < 5
    if ishandle(ii)
        count = count+1; 
        figure(ii)
        saveas(gcf, [save_dir filesep 'Seq_iter_' num2str(ii)], 'fig')
        pause(1)
        saveas(gcf, [save_dir filesep 'Seq_iter_' num2str(ii)], 'png')
    end
    end
end


end % end function. 