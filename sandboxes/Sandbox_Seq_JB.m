%% Sandbox_Seq_JB:
%
%
%
% EC 2021-03-10   initial version
%
%
%
%% get code and data paths
if exist('/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/Seq_JB', 'dir')
    addpath(genpath('/home/williamslab/Documents/Github/CEH2'))
    data_dir = '/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/Seq_JB';
    seq_dir = '/home/williamslab/Documents/Github/seqNMF'; 
else
    addpath(genpath('/home/ecarmichael/Documents/GitHub/CEH2'))
    data_dir = '/mnt/Data/Seq_JB';
    seq_dir = '/home/ecarmichael/Documents/GitHub/seqNMF';
end
cd(data_dir)

seq_data = load('M1269DMS1seqNMF.mat', 'seqNMF');
seq_data = seq_data.seqNMF;


addpath(seq_dir);

%% prepare data
place_cell_idx = seq_data.Place == 1;

cluster_id = seq_data.Kcluster;

Fs = seq_data.Cal.Dt;

data = seq_data.Cal.CellFreq';
% data = seq_data.Cal.Cell_Denoised'; 

movement_thresh = 2.5; % in cm/s
movement_idx = seq_data.Beh.TrackData.Speed >movement_thresh; % get times when the animal was moving.

pos(:,1) = seq_data.Beh.TrackData.X;
pos(:,2) = seq_data.Beh.TrackData.Y;

trial = seq_data.Beh.TrackData.Trial;

% start arm indices
SA_idx = ~ismissing(seq_data.Beh.TrackData.SA);

%reward indicies
R_idx = ~ismissing(seq_data.Beh.TrackData.EVENTS);

% get trial info as arrays
f_names = fieldnames(seq_data.Beh.GlobalActivity(1));
for iT = length(seq_data.Beh.GlobalActivity):-1:1
    for iF = 1:length(f_names)
        trial_times.(f_names{iF})(iT) = seq_data.Beh.GlobalActivity(iT).(f_names{iF});
    end
end

%  limit data to cluster groups

data_k1 = data(:,cluster_id == 1);
% [~, s_idx] = sort(data_k1,1);

data_k2 = data(:,cluster_id == 2);
data_k3 = data(:,cluster_id == 3);
data_k4 = data(:,cluster_id == 4);

%% get the SA per trial indicies

for iT = 1:max(unique(trial))
    all_trial_idx(iT,:) = SA_idx & trial(1:length(SA_idx)) == iT;
end


% Region_val = nan(1,length(SA_idx));
all_arms = []; count = 0; labels = [];
for iArm = unique(seq_data.Beh.TrackData.Region(~ismissing(seq_data.Beh.TrackData.Trial)));
    if contains(iArm(1), 'B_') || contains(iArm, 'G_')
        count = 1+ count;
        this_arm = zeros(1,length(seq_data.Beh.TrackData.Region));
        this_arm(contains(seq_data.Beh.TrackData.Region, iArm)) = 1;
        all_arms(count,:) = zeros(1,length(seq_data.Beh.TrackData.Region));
        all_arms(count,:)= logical(this_arm);
        labels{count} = strrep(iArm, '_', '');
    end
end
all_arms = logical(all_arms(:,1:length(SA_idx)));



%% try it with all the start arms togerther
block_idx = diff(trial(SA_idx));
block_idx = find(block_idx == 1);


data_in = data_k3(SA_idx,:)';
data_label = 'clust3';

%% sanity plots
% figure(100)
% subplot(1,3,3)
% hold on
% plot(pos(:,1), pos(:,2), '.k')
% c_ord = linspecer(size(all_arms,1));
%
% for iT = 1:size(all_arms,1)
%     if iT < 5
%         iC = 1;
%     elseif (4< iT)  && (iT < 8)
%         iC = 2;
%     elseif (8< iT)  && (iT < 12)
%         iC = 3;
%     elseif 12< iT
%         iC = 4;
%     end
%     plot(pos(all_arms(iT,:) & SA_idx,1), pos(all_arms(iT,:) & SA_idx,2), '.', 'color', c_ord(iT,:))
%     %     pause(1)
%     legend(labels)
%     vlines{iT} = 'w';
% end
% axis off
%
% %%
figure(101)
subplot(1,3,3)
hold on
plot(pos(:,1), pos(:,2), '.k')
c_ord = linspecer(size(all_trial_idx,1));
for iT = 1:size(all_trial_idx,1)
    plot(pos(all_trial_idx(iT,:),1), pos(all_trial_idx(iT,:),2), '.', 'color', c_ord(iT,:))
    % %     pause(1)
    these_labels{iT} = num2str(iT);
    %     vlines{iT} = 'w';
    
end
legend(['' these_labels])

axis off


subplot(1,3,1:2)
if length(unique(data_in)) < 3  % is this binary data.  May not be in logical format. 
    MS_Ca_Raster(data_in)
else
    MS_plot_ca_trace(data_in, [], 0.0001,1)
end
for iT = 1:size(all_trial_idx,1)-1
    h = vline(block_idx(iT), 'w', num2str(iT+1));
    h.Color  = c_ord(iT+1,:);
end

set(gcf, 'position', [212 846 1400 550])
%% try Seq on the merged SA blocks
% make a fake matrix of trial IDs for plotting.  This will be a matrix with
% trial ids across rows:
% 1 1 1 1 0 0 0 0 0 0 0 0 ...
% 0 0 0 0 2 2 2 2 0 0 0 0 ...
% 0 0 0 0 0 0 0 0 3 3 3 3 ...

trial_mat = trial(SA_idx);

trial_val_mat = zeros(max(unique(trial_mat)), size(trial_mat,2));
for ii = unique(trial_mat)
    trial_val_mat(ii, trial_mat == ii) = ii;
end

%% make a zero-padded version of the data to avoid sequences across start blocks.
% break data into trial blocks
this_idx = trial(SA_idx);
for iB = length(block_idx):-1:1
    trial_blocks{iB} = data_in(:,this_idx == iB);
    pos_blocks{iB} = trial_val_mat(:,this_idx == iB);
end


pad_data = []; pad_pos = []; % empty data and position matricies to fill.

pad = zeros(size(data_in,1), 20*Fs); % pad with 20s of zeros.
p_pad = zeros(size(trial_val_mat,1), 20*Fs);
pad_block = [];
for iB = 1:length(block_idx)
    if iB ==1
        pad_data = [trial_blocks{iB}, pad_data];
        pad_pos = [pos_blocks{iB}, pad_pos];
        pad_block(iB,1) = 1;
        pad_block(iB,2) = length(trial_blocks{iB});
    else
        
        pad_block(iB,1) = length([pad_data, pad])+1;
        pad_block(iB,2) = length([pad_data, pad, trial_blocks{iB}]);
        
        pad_data = [pad_data, pad, trial_blocks{iB}];
        pad_pos = [pad_pos, p_pad, pos_blocks{iB}];
    end
    %     pad_block(iB) = length(pad_pos) +1;
end



%%  prepare the data
%
% seq_data_in = data_in;
% seq_pos_in = trial_val_mat;

% using zer-pad
seq_data_in = pad_data;
seq_pos_in = pad_pos;


% break data into training set and test set
splitN = floor(size(seq_data_in,2)*.75);
splitS = floor(size(seq_pos_in,2)*.75);

% use the 12th trial as a split instead.
trial_split =  find(seq_pos_in(12,:) == 12);

% make the test/trian sets.
trainNEURAL = seq_data_in(:,1:trial_split(1)-1);
trainPOS = seq_pos_in(:,1:trial_split(1)-1);
testNEURAL = seq_data_in(:,(trial_split):end);
testPOS = seq_pos_in(:,(trial_split):end);

if ishandle(10)
    close(10)
end
h = figure(10);
h.Name = 'Seq input data splits';

subplot(2,4,1:3)

if length(unique(trainPOS)) < 3  % is this binary data.  May not be in logical format. 
MS_Ca_Raster(trainPOS, [], 4);
else
    MS_plot_ca_trace(trainPOS, [], 0.0001,1)
end
title('trial/position data (training)')


subplot(2,4,4)
if length(unique(testPOS)) < 3  % is this binary data.  May not be in logical format.
    MS_Ca_Raster(testPOS, [], 4);
else
    MS_plot_ca_trace(testPOS, [], 0.0001,1)
end
title('testing')
% axis off

subplot(2,4,5:7)
if length(unique(trainNEURAL)) < 3  % is this binary data.  May not be in logical format.
    MS_Ca_Raster(trainNEURAL, [], 4);
else
    MS_plot_ca_trace(trainNEURAL, [], 0.0001,1)
end
title('neural data (training)')

subplot(2,4,8)
if length(unique(testNEURAL)) < 3  % is this binary data.  May not be in logical format.
    MS_Ca_Raster(testNEURAL, [], 4);
else
    MS_plot_ca_trace(testNEURAL, [], 0.0001,1)
end% axis off
title('testing')

%% run SeqNMF across multiple time scales.


Ls = 2:12;

for iL = Ls
    % Set some parameters
    rng(235); % fixed rng seed for reproduceability
    X = trainNEURAL;
    K = 3;
    L = iL;
    Lneural = ceil(L*Fs);
    lambda = 0.001;% can be derived below. "Go just above the cross point"  this value is from prior testing.
    
    
    % all_seq.subject = subject;
    % all_seq.session = session;
    % all_seq.task = dir_parts{end}(end-3:end);
    %     all_seq.analysis_date = datetime;
    %     all_seq.data = data_in;
    % all_seq.pos = [];
    
    %     %% Procedure for choosing K
    %     tic
    %     Ws = {};
    %     Hs = {};
    %     numfits = 10; %number of fits to compare
    %     for k = 1:10
    %         display(sprintf('running seqNMF with K = %i',k))
    %         for ii = 1:numfits
    %             [Ws{ii,k},Hs{ii,k}] = seqNMF(X,'K',k, 'L', L,'lambda', lambda,'maxiter',100,'showplot',0);
    %             % note that max iter set low (30iter) for speed in demo (not recommended in practice)
    %         end
    %         inds = nchoosek(1:numfits,2);
    %         for i = 1:size(inds,1) % consider using parfor for larger numfits
    %             Diss(i,k) = helper.DISSX(Hs{inds(i,1),k},Ws{inds(i,1),k},Hs{inds(i,2),k},Ws{inds(i,2),k});
    %         end
    %
    %     end
    %     %% Plot Diss and choose K with the minimum average diss.
    %     figure,
    %     plot(1:10,Diss,'ko'), hold on
    %     h1 = plot(1:10,median(Diss,1),'k-','linewidth',2);
    %     % h2 = plot([3,3],[0,0.5],'r--');
    %     legend(h1, {'median Diss'})
    %     xlabel('K')
    %     ylabel('Diss')
    
    
    % %%  Find a good lambda
    % nLambdas = 30; % increase if you're patient
    % K = 2;
    % X = trainNEURAL;
    % lambdas = sort([logspace(-2,-4,nLambdas)], 'ascend');
    % loadings = [];
    % regularization = [];
    % cost = [];
    % [N,T] = size(X);
    %
    % for li = 1:length(lambdas)
    %     tic
    %     [W, H, ~,loadings(li,:),~]= seqNMF(X,'K',K,'L',Lneural,...
    %         'lambdaL1W', .1, 'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0);
    %     [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X,W,H);
    %     display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
    %     toc
    % end
    % %% plot Lambda cost
    %
    % windowSize = 3;
    % b = (1/windowSize)*ones(1,windowSize);
    % a = 1;
    % Rs = filtfilt(b,a,regularization);
    % minRs = prctile(regularization,10); maxRs= prctile(regularization,90);
    % Rs = (Rs-minRs)/(maxRs-minRs);
    % R = (regularization-minRs)/(maxRs-minRs);
    % Cs = filtfilt(b,a,cost);
    % minCs =  prctile(cost,10); maxCs =  prctile(cost,90);
    % Cs = (Cs -minCs)/(maxCs-minCs);
    % C = (cost -minCs)/(maxCs-minCs);
    %
    % figure; hold on
    % plot(lambdas,Rs, 'b')
    % plot(lambdas,Cs,'r')
    % scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
    % scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
    % xlabel('Lambda'); ylabel('Cost (au)')
    % set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
    % set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
    % set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
    % shg
    
    
    
    %% choose lambda=.005; run multiple times, see number of sig factors
    %     loadings = [];
    %     pvals = [];
    %     is_significant = [];
    %     X = trainNEURAL;
    %     nIter = 20; % increase if patient
    % %     K = 3;
    %
    %     display(['Running seqNMF multiple times for lambda=' num2str(lambda)])
    %
    %     for iteri = nIter:-1:1
    %         display(['seqNMF run ' num2str(iteri) '/' num2str(nIter)])
    %         [W, H, ~,loadings(iteri,:),power]= seqNMF(X,'K',K,'L',Lneural,...
    %             'lambdaL1W', .1, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0);
    %         p = .05;
    %         [pvals(iteri,:),is_significant(iteri,:)] = test_significance(testNEURAL,W,p);
    %         if sum(is_significant) ==0
    %             continue
    %         end
    %         W = W(:,is_significant(iteri,:)==1,:);
    %         H = H(is_significant(iteri,:)==1,:);
    %         [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
    %         indSort = hybrid(:,3);
    %         tstart = 300;
    %         clf; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 0,trainPOS(:,floor(tstart*Fs/Fs):end))
    %         display(['Sig factors found: ' num2str(sum(is_significant(iteri,:)))])
    %
    %     end
    %     figure; hold on
    %     h = histogram(sum(is_significant,2), 'edgecolor', 'w', 'facecolor', .7*[1 1 1]);
    %     h.BinCounts = h.BinCounts/sum(h.BinCounts)*100;
    %     xlim([0 10]);
    %     xlabel('# significant factors')
    %     ylabel('% seqNMF runs')
    
    
    %% PARTS choose lambda=.005; run multiple times, see number of sig factors
    loadings = [];
    pvals = [];
    is_significant = [];
    X = trainNEURAL;
    nIter = 30; % increase if patient
    
    display(['Running seqNMF multiple times for lambda=' num2str(lambda)])
    lambdaOrthoH = 0;
    lambdaOrthoW = 1; % favor parts-based (these can take any value, don't need to be zero and one)
    for iteri = nIter:-1:1
        tic
        display(['seqNMF run Parts ' num2str(iteri) '/' num2str(nIter)])
        [W, H, ~,loadings(iteri,:),power]= seqNMF(X,'K',K,'L',Lneural,...
            'lambdaOrthoH', lambdaOrthoH, 'lambdaOrthoW', lambdaOrthoW, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0);
        p = .05;
        [pvals(iteri,:),is_significant(iteri,:)] = test_significance(testNEURAL,W,p);
        if sum(is_significant) ==0
            continue
        end
        W = W(:,is_significant(iteri,:)==1,:);
        H = H(is_significant(iteri,:)==1,:);
        if isempty(W)
            continue
        end
        all_W{iteri} = W;
        all_H{iteri} = H;
        [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
        indSort = hybrid(:,3);
        tstart = 1;
        figure(iteri)
        WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), 0,trainPOS(:,floor(tstart*Fs/Fs):end))
        display(['Sig factors found: ' num2str(sum(is_significant(iteri,:)))])
        toc
        
            close all
    end
    
    
    %% saver evertyhing
    Seq_out{iL}.K = K;
    Seq_out{iL}.L = L;
    Seq_out{iL}.W = all_W;
    Seq_out{iL}.H = all_H;
    Seq_out{iL}.is_significant  =is_significant;
    Seq_out{iL}.lambda = lambda;
    Seq_out{iL}.lambdaOrthoH = lambdaOrthoH;
    Seq_out{iL}.lambdaOrthoW = lambdaOrthoW;
    save([ data_dir filesep data_label '_Seq_out_pad.mat'], 'Seq_out');
    
    %% plot the raw data with a significant sequence
    
    if sum(Seq_out{iL}.is_significant,'all')
        
        % histogram
        figure(200+iL); hold on
        h = histogram(sum(Seq_out{iL}.is_significant,2), 'edgecolor', 'w', 'facecolor', .7*[1 1 1]);
        h.BinCounts = h.BinCounts/sum(h.BinCounts)*100;
        xlim([0 10]);
        xlabel('# significant factors')
        ylabel('% seqNMF runs')
        set(gcf, 'position', [212 846 1400 550])
        saveas(gcf, [data_dir filesep data_label '_K_hist_' num2str(Seq_out{iL}.K) '_L' num2str(iL) '.fig'])
        saveas(gcf, [data_dir filesep data_label '_K_hist_' num2str(Seq_out{iL}.K) '_L' num2str(iL) '.png'])
        close(200+iL)
        
        
        if max(sum(Seq_out{iL}.is_significant,2)) > 1
            [~,sig_iter] = max(sum(Seq_out{iL}.is_significant,2));
        else
            sig_iter = find(sum(Seq_out{iL}.is_significant,2));
        end
        
        this_sig_iter = sig_iter(1);
        s_ord = linspecer(max(sum(Seq_out{iL}.is_significant,2)));
        
        [~, ~, ~, hybrid] = helper.ClusterByFactor(Seq_out{iL}.W{this_sig_iter}(:,:,:),1);
        indSort = hybrid(:,3);
        
        % make a figure
        figure(110+iL)
        % info
        subplot(5,8,1)
        text(0, 0.75, ['L = ' num2str(iL)]);
        text(0, 0.25, ['iter = ' num2str(this_sig_iter)]);
        axis off
        
        % plot the trials
        subplot(5,8,[8, 16, 24, 32, 40]);
        hold on
        plot(pos(:,1), pos(:,2), '.k')
        c_ord = linspecer(size(all_trial_idx,1));
        for iT = 1:size(all_trial_idx,1)
            plot(pos(all_trial_idx(iT,:),1), pos(all_trial_idx(iT,:),2), '.', 'color', c_ord(iT,:))
            % %     pause(1)
            these_labels{iT} = num2str(iT);
            %     vlines{iT} = 'w';
        end
        %         legend(['' these_labels], 'location', 'southoutside', 'orientation', 'horizontal')
        axis off
        
        % plot the H factors
        if size(Seq_out{iL}.H{this_sig_iter},1) == 1
            ax1= subplot(5,8,2:7);
        elseif size(Seq_out{iL}.H{this_sig_iter},1) == 2
            ax1=  subplot(5,8,3:7);
        elseif size(Seq_out{iL}.H{this_sig_iter},1) == 3
            ax1=  subplot(5,8,4:7);
        end
        for ii = 1:size(Seq_out{iL}.H{this_sig_iter},1) % loop sig factors
            hold on
            plot((0:size(Seq_out{iL}.H{this_sig_iter},2)-1)/Fs, Seq_out{iL}.H{this_sig_iter}(ii,:)+ii*.5, 'color', s_ord(ii,:))
        end
        set(gca, 'ytick', 1:size(all_trial_idx,1))
        
        
        % plot W factors
%         if size(Seq_out{iL}.H{this_sig_iter},1) > 1 % only works for 2 factors b/c lazy.
            for ii = 1:size(Seq_out{iL}.H{this_sig_iter},1)
                if ii == 1
                    subplot(5,8,[9 17 25 33])
                elseif ii == 2
                    subplot(5,8,[10 18 26 34])
                    ylabel([]);
                elseif ii == 3
                    subplot(5,8,[11 19 27 35])
                    ylabel([]);
                end
                this_seq = squeeze(Seq_out{iL}.W{this_sig_iter}(indSort,ii,:));
                this_seq(this_seq > max(this_seq)*0.75) = 1;
                this_seq(this_seq ~=1) = 0;
                MS_Ca_Raster(this_seq,(0:size(this_seq,2)-1)/Fs,2)
                %                 imagesc((0:size(this_seq,2)-1)/Fs, 1:size(this_seq,1), this_seq.*(1:size(this_seq,1))');
                xlabel('time (s)');
                title(['Factor ' num2str(ii) '(top 25%)']);
                %                 colormap(linspecer(size(this_seq,1)+1));
                %                 caxis([0 size(this_seq,1)]);
                %                 j_ord =jet(100);
                %                 set(gca, 'color', j_ord(1,:));
                set(gca, 'XColor', s_ord(1,:),'YColor', s_ord(1,:));
                if ii >1; set(gca, 'yticklabel', [], 'XColor', s_ord(ii,:),'YColor', s_ord(ii,:));  end
            end
%         else
%             subplot(5,8,[9 17 25 33])
%              this_seq = squeeze(Seq_out{iL}.W{this_sig_iter}(indSort,:,:));
%                 this_seq(this_seq > max(this_seq)*0.75) = 1;
%                 this_seq(this_seq ~=1) = 0;
% %             this_seq = squeeze(Seq_out{iL}.W{this_sig_iter}(indSort,:,:));
%                             MS_Ca_Raster(this_seq,(0:size(this_seq,2)-1)/Fs,2)
% %             imagesc((0:size(this_seq,2)-1)/Fs, 1:size(this_seq,1), this_seq.*(1:size(this_seq,2)));
%             xlabel('time (s)')
%             set(gca, 'XColor', s_ord(1,:),'YColor', s_ord(1,:))
% %             colormap(jet);
%         end
        %         ylim([0 size(all_H{this_sig_iter},1)+1])
        
        
        % plot the raw data sorted based on these factor
        if size(Seq_out{iL}.H{this_sig_iter},1) == 1
            ax2 = subplot(5,8,[10:15 18:23 26:31 34:39]);
        elseif size(Seq_out{iL}.H{this_sig_iter},1) == 2
            ax2 = subplot(5,8,[11:15 19:23 27:31 35:39]);
        elseif size(Seq_out{iL}.H{this_sig_iter},1) == 3
            ax2 = subplot(5,8,[12:15 20:23 28:31 36:39]);
        end
        
        
        if length(unique(seq_data_in)) < 3  % is this binary data.  May not be in logical format. 
            MS_Ca_Raster(seq_data_in(indSort,:), (0:size(seq_data_in,2)-1)/Fs, 2)
        else
            MS_plot_ca_trace(seq_data_in(indSort,:),(0:size(seq_data_in,2)-1)/Fs, 0.0001, 1); 
        end
        %         imagesc((0:size(seq_data_in,2)-1)/Fs, 1:size(seq_data_in,1),seq_data_in(indSort,:))
        for iT = 1:size(all_trial_idx,1)-1
            rectangle('position', [pad_block(iT,1)/Fs, 0,(pad_block(iT,2)-pad_block(iT,1))/Fs , 2], 'facecolor', c_ord(iT+1,:))
            
            %             h = vline(block_idx(iT)/Fs, 'w', num2str(iT));
            %             h.Color  = c_ord(iT+1,:);
        end
        
        linkaxes([ax1 ax2], 'x')
        
        set(gcf, 'position', [212 846 1400 550])
        saveas(gcf, [data_dir filesep data_label '_Seq_parts_K' num2str(Seq_out{iL}.K) '_L' num2str(iL) '.fig'])
        saveas(gcf, [data_dir filesep data_label '_Seq_parts_K' num2str(Seq_out{iL}.K) '_L' num2str(iL) '.png'])
        close(110+iL)
        
    end
    
    
    
    close all
    
    
end

%%
for iL = Ls
    %     if sum(Seq_out{iL}.is_significant, 'all')
    %     fprintf('Sig factors (%i/%i; 1 = %i, 2 = %i, 3 = %i) found in L = %i\n', sum(Seq_out{iL}.is_significant, 'all'),nIter, iL)
    fprintf('Significant Seqs %d/%d inters (K1: %d, K2: %d, K3: %d) found using K: %d L: %ds\n',sum((sum(Seq_out{iL}.is_significant,2)~=0), 'all'),nIter, ...
        sum((sum(Seq_out{iL}.is_significant,2)==1),'all'), sum((sum(Seq_out{iL}.is_significant,2)==2),'all'), sum((sum(Seq_out{iL}.is_significant,2)==3),'all'),...
        Seq_out{iL}.K, iL);
    %     end
end

%% temp plotting stuff (remove)

hold on
for iT = 1:size(all_trial_idx,1)-1
    rectangle('position', [pad_block(iT,1)/Fs, 0,(pad_block(iT,2)-pad_block(iT,1))/Fs , 2], 'facecolor', c_ord(iT+1,:))
    
    % h = vline(pad_block(iT)/Fs, 'w', num2str(iT));
    % h.Color  = c_ord(iT+1,:);
    %
    % h = vline(pad_block(iT)+(20*Fs)/Fs, 'w', num2str(iT));
    % h.Color  = c_ord(iT+1,:);
end


%% rank trials based on H values

iL = 2;
if max(sum(Seq_out{iL}.is_significant,2)) > 1
    [~,sig_iter] = max(sum(Seq_out{iL}.is_significant,2));
else
    sig_iter = find(sum(Seq_out{iL}.is_significant,2));
end

this_sig_iter = sig_iter(1);

H = Seq_out{iL}.H(this_sig_iter);

s_ord = linspecer(max(sum(Seq_out{iL}.is_significant,2)));

[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(Seq_out{iL}.W{this_sig_iter}(:,:,:),1);
indSort = hybrid(:,3);


% split the H values into trials
%         for iT = unique(trial_blocks)



