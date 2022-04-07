function CC_SeqNMF_Sleep(data_in,Fs, k, l, lambda, save_dir )



if nargin < 2
    Fs = 33;
    k = 2;
    l = 10;
    lambda = 0.005;
    save_dir = [cd filesep 'Seq_out_sleep_L' strrep(num2str(l), '.', 'p')];
elseif nargin < 3
    k = 2;
    l = 10;
    lambda = 0.005;
    save_dir = [cd filesep 'Seq_out_sleep_L' strrep(num2str(l), '.', 'p')];
elseif nargin < 4
    lambda = 0.005;
    l = 10;
    save_dir = [cd filesep 'Seq_out_sleep_L' strrep(num2str(l), '.', 'p')];
elseif nargin < 5
    lambda = 0.005;
    save_dir = [cd filesep 'Seq_out_sleep_L' strrep(num2str(l), '.', 'p')];
elseif nargin < 6
    save_dir = [cd filesep 'Seq_out_sleep_L' strrep(num2str(l), '.', 'p')];
end


%% Fit with seqNMF: most of this is straight out of Mackevicius et al. 2019 https://elifesciences.org/articles/38471

fprintf('Running all laps\n')
this_data = data_in;


% break data into training set and test set
splitN = floor(size(this_data,2)*.75); 
 
% make the test/trian sets.  
trainNEURAL = this_data(:,1:splitN); 
testNEURAL = this_data(:,(splitN+1):end); 


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
        figure(iteri); clf; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 1)
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
Seq_out_sleep.k = k;
Seq_out_sleep.l = l;
Seq_out_sleep.lambda = opt_lambda; 
Seq_out_sleep.lambdaOrthoH = lambdaOrthoH;
Seq_out_sleep.lambdaOrthoW = lambdaOrthoW; 
Seq_out_sleep.this_data = this_data; 
Seq_out_sleep.trainNEURAL = trainNEURAL; 
Seq_out_sleep.testNEURAL = testNEURAL; 
Seq_out_sleep.W = Ws;
Seq_out_sleep.H = Hs;
Seq_out_sleep.hybrid = hybrid; 
Seq_out_sleep.loadings = loadings;
Seq_out_sleep.is_significant = is_significant;
Seq_out_sleep.pvals = pvals; 

%% save come ouputs
mkdir(save_dir)
save([save_dir filesep 'Seq_out_sleep'], 'Seq_out_sleep',  '-v7.3'); 

figure(1010)
saveas(gcf, [save_dir filesep 'Sig_seqs_sleep_bar'], 'fig')
pause(1)
saveas(gcf, [save_dir filesep 'Sig_seqs_sleep_bar'], 'png')

count = 0; 
for ii = 1:10
    if count < 5
    if ishandle(ii)
        count = count+1; 
        figure(ii)
        saveas(gcf, [save_dir filesep 'Seq_sleep_iter_' num2str(ii)], 'fig')
        pause(1)
        saveas(gcf, [save_dir filesep 'Seq_sleep_iter_' num2str(ii)], 'png')
    end
    end
end


end % end function. 