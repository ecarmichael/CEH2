function MS_Run_REM_SeqNMF(cfg_in, data_dir, pro_dir, REM_dir)
%% MS_Run_SeqNMF: runs SeqNMF on a data session and the corresponding REM session. 
%
%  INTENDED TO BE CALLED IN Master_SeqNMF_JC
%
%    Inputs: 
%    - cfg_in :   [struct]  contains configuration paramaters which will
%    override the defaults
%
%    - data_dir:   [string]   path to data (ms.mat and behav.mat)
%
%    - pro_dir:  [string] path to the processed data (place cell
%    classification/anxiety...
%
%    -  REM_dir:  [string] path to REM data corresponding to data_dir. 
%
%
%    Outputs: 
%   
%
%
%
%
% EC 2021-02-25   initial version 
%
%
%
%% initialize



global PARAMS



% set up configuration
cfg_def = [];
cfg_def.place = 1; % run using only place cells
cfg_def.anx = 0; % run using only anxiety cells.
cfg_def.select_lambda = 0; 


cfg_def.score = 0; % run the sequence scoring.  Warning very very slow
cfg_def.score_nShuf = 500; 


% REM
cfg_def.REM = 1; % run the REM data



cfg = ProcessConfig2(cfg_def, cfg_in); 

if nargin < 4
    fprintf('No REM data specified!  only analyzing wake data\n')
    cfg.REM = 0; % flag to skip REM later on. 
end


if cfg.place && ~cfg.anx
    type = 'Place';
    fprintf('Running for Place cells only\n');
elseif ~cfg.place && cfg.anx
        type = 'anx';
    fprintf('Running for Anxiety cells only\n');
elseif ~cfg.place && ~cfg.anx
        type = 'all';
    fprintf('Running for all cells\n');
elseif cfg.place && cfg.anx
    error('Can''t run place and anxiety in one go.  Run for each separately.')
end
%% load a nice session

dir_parts = strsplit(data_dir, filesep);
task = strsplit(dir_parts{end}, '_');
task = task{end}; 
if contains(task, 'HATD6')
    task = 'HATDSwitch';
end
subject = dir_parts{end-2};
session = dir_parts{end-1};


if cfg.REM
% REM file dir
this_rem_dir = [REM_dir filesep upper(subject) filesep session];
end

% load the main data
cd(data_dir);

warning off  % otherwise you get all the videoencoder warnings.
load('ms.mat');
% clearvars -except ms this_sess PARAMS  % some ms files saved all the other variables as well.  not sure why.
load('behav.mat');
if exist('AllSpatialFiringData.mat', 'file')
    load('AllSpatialFiringData.mat')
    fprintf('<strong>No AllSpatialFiringData.mat found.  Can''t distinguish place cells</strong>\n')
end

if ~isfield(ms, 'Binary')
    fprintf('No binary field found in ms, extracting...\n')
    ms = MS_msExtractBinary_detrendTraces(ms, 2);
end

% use another place cell classification file
cd(pro_dir)

if contains(data_dir, 'HAT')
    if exist([pro_dir  '3.Anxiety_SafetyCell' filesep subject filesep 'Anxiety_output' filesep 'STD1'], 'dir')
        load([pro_dir '3.Anxiety_SafetyCell' filesep subject filesep 'Anxiety_output' filesep 'STD2' filesep 'Anxiety_Output_std2.mat'])
    else
        load([pro_dir  '3.Anxiety_SafetyCell' filesep subject filesep 'Anxiety_output'  filesep 'Anxiety_Output_std2.mat'])
    end
    % get place cells as well. 
    if contains(task, 'HATD')% workaround for naming of HAT + D
        load([pro_dir filesep '4.PlaceCell' filesep subject filesep task filesep 'spatial_analysis.mat'])
        load([pro_dir filesep '4.PlaceCell' filesep subject filesep task filesep 'SA.mat'])
    else
        load([pro_dir filesep '4.PlaceCell' filesep subject filesep strrep(task, 'HAT', 'HATD') filesep 'spatial_analysis.mat'])
        load([pro_dir filesep '4.PlaceCell' filesep subject filesep strrep(task, 'HAT', 'HATD') filesep 'SA.mat'])
    end
else
    load("spatial_analysis.mat")
    load("spatial_analysis_classif.mat")
    fprintf('Using Processed spatial info file ''spatial_analyses_classif.mat''\n')
end
%% plug it into SeqNMF

addpath(PARAMS.code_seqnmf_dir)
Fs =mode(diff(ms.time));

% data_in = ms.FiltTraces';
data_in = ms.Binary';

pos(:,1) = interp1(behav.time,behav.position(:,1),ms.time);
pos(:,2) = interp1(behav.time,behav.position(:,2),ms.time);
velocity = interp1(behav.time, behav.speed, ms.time);

% remove inactive cells
% total_act = sum(data_in,2);

% keep_idx = total_act>30; % exclude cells that are active for less than 1s total

% data_in = data_in(keep_idx,:);

% sort data based on place field location

% centroids = [];
% for iC = length(spatial_analysis.bin):-1:1
%     if iscell(spatial_analysis.bin{iC,1}.PlaceFieldCentroid)
%         temp = cell2mat(spatial_analysis.bin{iC,1}.PlaceFieldCentroid);
%     else
%         temp = spatial_analysis.bin{iC,1}.PlaceFieldCentroid;
%     end
%     centroids(iC) = temp(1);
%
%     % is it a place cell?
%     Place_Cell_D1(iC) = spatial_analysis.bin{iC,1}.IsPlaceCell;
%
%     place_map_1d(iC,:) = mean(spatial_analysis.bin{iC,1}.PlaceField,1)/max(mean(spatial_analysis.bin{iC,1}.PlaceField,1)); % get the 1d place field and normalize.
% end

% limit to predefined anxiety cells.
place_cell_idx  = SA.WholePlaceCell; 
if contains(data_dir, 'HAT'); anx_cell_idx = (Anxiety_Output.AnxietyCell_day_switch); end


% for labeling saved files later. 
if cfg.place
    keep_idx = sort(place_cell_idx); 
    centroids = NaN*zeros(1,size(data_in,1)); 
    place_maps = NaN*zeros(size(data_in,1),20);  % assumes 20bins.  fit for future. 
    for iC = length(centroids):-1:1
        if ismember(iC, keep_idx)
            for iT = 3:-1:1
                if spatial_analysis.bin{iC,iT}.IsPlaceCell
                    centroids(iC) = spatial_analysis.bin{iC,iT}.PlaceFieldCentroid{1}(1);
                    place_maps(iC,:) = spatial_analysis.bin{iC, iT}.PlaceField'; 
                end
            end
        end
    end
    place_maps(isnan(centroids),:) = [];
    centroids(isnan(centroids)) = []; 
elseif cfg.anx
    keep_idx = anx_cell_idx;
    fprintf('%0.0f%% of anxiety cells are also place cells\n', 100*sum(ismember(anx_cell_idx, place_cell_idx))/length(anx_cell_idx))
else
    keep_idx = 1:size(ms.Binary,2); % just use all the cells
end

    fprintf('Running analyses on %s cells %d/%d = %0.0f%% ...\n', type, length(keep_idx), size(data_in,1),(length(keep_idx)/size(data_in,1))*100)


keep_idx(isnan(keep_idx)) = [];
keep_idx = sort(keep_idx);
data_in = data_in(keep_idx,:);

% sort the centroids (1D only)
%     centroids = centroids(Place_Cell_D1);
%
% [~,sort_idx] = sort(centroids);
%
% place_map_1d = place_map_1d(Place_Cell_D1,:);
% place_map_1d = place_map_1d(sort_idx,:);
%
% imagesc(place_map_1d);
% axis xy

% data_in = data_in(sort_idx,:);

% %% run the normalizing steps from SeqNMFE paper: Since seqNMF would
% % prioritize the neurons with the most power, we renormalized by dividing
% % the signal from each neuron by the sum of the maximum value of that row
% % and the 95th percentile of the signal across all neurons. In this way,
% % neurons with larger peaks were given some priority, but not much more
% % than that of neurons with weaker signals.
%
% p_95 = prctile(data_in, 95, 'all');
% for iC = size(data_in, 1):-1:1
%    data_in(iC,:) =  (data_in(iC,:) / sum(data_in(iC,:)))/p_95;
%    data_in(iC,:) = data_in(iC,:) + abs(min(data_in(iC,:)));
%
% end
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

% figure(111)
% hold on
% % plot(ms.time, pos, 'color', [.2 .2 .2])
% plot(L_laps_time/1000, nan(length(L_laps_pos),1), 'b')
% plot(R_laps_time/1000, nan(length(R_laps_pos),1), 'r')
% 
% plot(L_laps_time/1000, L_laps_pos, '.b')
% plot(R_laps_time/1000, R_laps_pos, '.r')
% xlabel('time (s)')
% ylabel('position')
% legend({'Left Laps', 'Right Laps'})
% title('Position (L/R laps)')

% split the velocity values
L_velo = []; R_velo = [];
for iL = 1:length(L_lap_start_idx)
    L_velo = [L_velo; velocity(L_lap_start_idx(iL):L_lap_end_idx(iL))];% append the velocity values for each R lap
end
for iL = 1:length(R_lap_start_idx)
    R_velo = [R_velo; velocity(R_lap_start_idx(iL):R_lap_end_idx(iL))]; % append the velocity values for each R lap
end

% convert trajectory to 3d mat as place x time x logical (was the animal here) for plotting
% in the SeqNMFe WHPlot function
p_range = floor(min(L_laps_pos(:,1),[], 'all')):0.5:ceil(max(L_laps_pos,[],'all')); % make a range of positions which the middle value here being the precision (default 0.5cm)
L_pos_mat = [];
for iP = length(L_laps_pos):-1:1
    pos_idx= nearest_idx3(L_laps_pos(iP), p_range);
    pos_val = zeros(size(p_range));   % hack to have zeros in all the other positions.
    pos_val(pos_idx) = 1;
    %     pos_val(pos_idx) = L_velo(iP);
    L_pos_mat(iP,:) = smooth(pos_val,10); % fill in the y axis with zeros for unoccupied and 1 for occupied.
end

p_range = floor(min(R_laps_pos(:,1),[], 'all')):0.5:ceil(max(R_laps_pos,[],'all')); % make a range of positions which the middle value here being the precision (default 0.5cm)
R_pos_mat = [];
for iP = length(R_laps_pos):-1:1
    pos_idx= nearest_idx3(R_laps_pos(iP), p_range);
    pos_val = zeros(size(p_range));   % hack to have zeros in all the other positions.
    %     pos_val(pos_idx) = R_velo(iP);
    pos_val(pos_idx) = 1;
    R_pos_mat(iP,:) = smooth(pos_val,10); % fill in the y axis with zeros for unoccupied and 1 for occupied.
end

% make a position x time matrix for plotting
p_range - floor(min(pos, [], 'all')) : 0.5: ceil(max(pos,[],'all'));
pos_mat = [];
for iP = length(pos):-1:1
    if isnan(pos(iP))
        pos_mat(iP,:) = zeros(size(p_range));
    else
        pos_idx= nearest_idx3(pos(iP), p_range);
        pos_val = zeros(size(p_range));   % hack to have zeros in all the other positions.
        %     pos_val(pos_idx) = R_velo(iP);
        if pos_idx >= length(p_range)-3
            pos_val(pos_idx-3:end) = iP;
        elseif pos_idx <= 3
            pos_val(1:pos_idx+3) = iP;
        else
            pos_val(pos_idx-3:pos_idx+3) = iP;
        end
        %         pos_val(pos_idx) = iP;
        pos_mat(iP,:) = pos_val; % fill in the y axis with zeros for unoccupied and 1 for occupied.
    end
end
pos_mat= pos_mat';
% figure(112);
% plot_pos_mat = pos_mat;
% plot_pos_mat(plot_pos_mat == 0) =NaN;
% % plot_pos_mat(~isnan(plot_pos_mat)) = 1;
% nan_imagesc_ec(plot_pos_mat);
% axis xy
% x_ticks = get(gca, 'xtick');
% set(gca, 'xticklabels', x_ticks/Fs)
% set(gca, 'color', 'w')
% colormap('jet')
% title('position matrix')


if cfg.place
    figure(10)
    subplot(5,4,1)
    parts = strsplit(session, 'P');
    text(0, .8, (strrep(parts{1}, '_', ' ')))
    text(0, .4, ['P' (strrep(parts{2}, '_', ' '))])
    axis off;
    
    ax1(1) =  subplot(5,4,2:4);
    temp_pos = pos_mat;
    temp_pos(temp_pos == 0) = NaN;
    h = nan_imagesc_ec(temp_pos);
    h.XData = (1:length(pos_mat)+1)/Fs;
    axis xy
    x_ticks = get(gca, 'xtick');
    set(gca, 'xticklabels', x_ticks/Fs)
    set(gca, 'color', 'w')
    %     imagesc((1:size(this_pos, 2))/Fs,1:size(this_pos,1),this_pos.*([1:size(this_pos,2)]+500))
    set(gca, 'xtick', [], 'ytick', []);
    title('position')
    
    ax2(1) = subplot(5,4,[5 9 13 17]);
    [~, idx] = sort(centroids);
    this_seq = place_maps(idx,:);
    this_seq(this_seq > 0) = 1;
    imagesc((1:size(this_seq, 2))/Fs,1:size(this_seq, 1), this_seq.*([1:size(this_seq,1)]+(floor(size(this_seq,1)/2)))');
    ylabel('cell ID (sorted)')
    xlabel('seq time (s)');
    
    ax1(2) = subplot(5,4,[6:8 10:12 14:16 18:20]);
    imagesc((1:size(data_in, 2))/Fs,1:size(data_in, 1), data_in(idx,:).*([1:size(data_in,1)]+(floor(size(data_in,1)/2)))');
    set(gca, 'ytick', []);
    xlabel('concatenated time (s)')
    cfg_plot.ft_size = 12;
    SetFigure(cfg_plot, gcf);
    set(gcf, 'position', [680 500 1127 480])
    linkaxes(ax1, 'x')
    
    colormap('jet');
    
    colormap(ax2(1), 'parula');
    colormap(ax1(2), 'parula');
    
    saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Raw_sorted' '_' type])
    saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Raw_sorted' '_' type '.png'])
    
    %     cmap = linspecer(size(data_in, 1));
    % %         cmap(1,:) = 0;
    % colormap(gca, cmap);
    
    % get some zoomed in plots. 
    for ii = floor(max((1:size(data_in, 2))/Fs)/6): floor(max((1:size(data_in, 2))/Fs)/3):max((1:size(data_in, 2))/Fs)
        xlim([ii-floor(max((1:size(data_in, 2))/Fs)/6) ii])

        
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Raw_sorted_zoom' num2str(ii-floor(max((1:size(data_in, 2))/Fs)/6)) '_' type])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Raw_sorted_zoom' num2str(ii-floor(max((1:size(data_in, 2))/Fs)/6)) '_' type '.png'])
        %
    end
    xlim([0 max((1:size(data_in, 2))/Fs)])
end
%% Fit with seqNMF: most of this is straight out of Mackevicius et al. 2019 <https://elifesciences.org/articles/38471 https://elifesciences.org/articles/38471>

% this_data = data_in; % which data to use.
% this_data = R_data;
this_data = data_in(:,movement_idx);
% this_pos = R_pos_mat';
this_pos = pos_mat(:,movement_idx);

PosFs = Fs; % b/c we already interpolated the position data to the ms.time.

% break data into training set and test set
splitN = floor(size(this_data,2)*.75);
splitS = floor(size(this_pos,2)*.75);

% make the test/trian sets.
trainNEURAL = this_data(:,1:splitN);
trainPOS = this_pos(:,1:splitS);
testNEURAL = this_data(:,(splitN+1):end);
testPOS = this_pos(:,(splitS+1):end);

% Set some parameters
rng(235); % fixed rng seed for reproduceability
X = trainNEURAL;
K = 3;
L = 10;
Lneural = ceil(L*Fs);
    lambda = 0.001;% from the above.  "Go just above the cross point"  this value is from prior testing.  Can be recomputed using cfg.select_lambda = 1; 


all_seq.subject = subject;
all_seq.session = session;
all_seq.task = dir_parts{end}(end-3:end);
all_seq.analysis_date = datetime;
all_seq.data = this_data;
all_seq.pos = this_pos;

%% Procedure for choosing lambda

if cfg.select_lambda
nLambdas = 30; % increase if you're patient
K = 3;
X = trainNEURAL;
lambdas = sort([logspace(-2,-4,nLambdas)], 'ascend');
loadings = [];
regularization = [];
cost = [];
    [N,T] = size(X);

parfor li = 1:length(lambdas)
    [W, H, ~,loadings(li,:),~]= seqNMF(X,'K',K,'L',Lneural,...
        'lambdaL1W', .1, 'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0);
    [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X,W,H);
    display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
end
%% plot Lambda cost

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

figure; hold on
plot(lambdas,Rs, 'b')
plot(lambdas,Cs,'r')
scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
xlabel('Lambda'); ylabel('Cost (au)')
set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])
shg

% get the lambda point after the cross point.
lambda_diff = Rs - Cs;
idx = find(lambda_diff < 0);

lambda = lambdas(idx(1));
end
%% Events based version

fprintf('\n*************** Events based ***************\n')
type = 'events';
% event based
lambdaOrthoH = .1; % favor events-based (these can take any value, don't need to be zero and one)
lambdaOrthoW = 0;
fprintf('Running seqNMF...K = %d  L = %d sec\n', K, L)
figure; 
[W_e, H_e, ~,~,~]= seqNMF(X,'K',K,'L',Lneural,...
    'lambdaL1W', .1, 'lambda', lambda, 'maxiter', 100, 'showPlot', 1,...
    'lambdaOrthoH', lambdaOrthoH, 'lambdaOrthoW', lambdaOrthoW);
% test if any are significant
p = .05; % desired p value for factors

disp('Testing significance of factors on held-out data')
[~,is_significant] = test_significance(testNEURAL,W_e,p);

W_e = W_e(:,is_significant,:);
H_e = H_e(is_significant,:);
fprintf('Events-based: Found %d/%d significant factors\n', sum(is_significant), length(is_significant))

% Plot events-based results
if isempty(H_e)
    fprintf('<strong>No significant sequences detected.  Skipping...</strong>\n')
else
    [~, ~, ~, hybrid] = helper.ClusterByFactor(W_e(:,:,:),1);
    indSort = hybrid(:,3);
    tstart = 1; % plot data starting at this timebin
    figure;
    WHPlot(W_e(indSort,:,:),H_e(:,tstart:end), X(indSort,tstart:end), ...
        0,trainPOS(:,floor(tstart*Fs/PosFs):end))
    title('Significant seqNMF factors, with raw data: Event-based')
    
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_fac_events' num2str(K) 'L' num2str(L) '_' type])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_fac_events' num2str(K) 'L' num2str(L) '_' type '.png'])
    
    figure;
    WHPlot(W_e(indSort,:,:),H_e(:,tstart:end), ...
        helper.reconstruct(W_e(indSort,:,:),H_e(:,tstart:end)),...
        0,trainPOS(:,floor(tstart*Fs/PosFs):end))
    title('SeqNMF reconstruction: Event-based')
    
     saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_recon_events' num2str(K) 'L' num2str(L) '_' type])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_recon_events' num2str(K) 'L' num2str(L) '_' type '.png'])
        
    addpath(genpath([PARAMS.code_seqnmf_dir filesep 'misc_elm']));
    figure; HTriggeredSpec(H_e,trainPOS,Fs,Fs,ceil(L*Fs));
    title('Seq triggerred position')
    
    figure; HTriggeredRaster(H_e,trainNEURAL(indSort,:),ceil(L*Fs));
    title('Seq triggered raster')
    
    figure(304)
    for iSeq = 1:size(W_e,2)
        subplot(1, size(W_e,2),iSeq)
        % offset = 10^-1; % offset for each cell
        this_seq = squeeze(W_e(indSort,iSeq,:));
        imagesc(1:size(this_seq,1)/Fs,1:size(this_seq,2), this_seq);
        ylabel('sorted cell ID')
        xlabel('seq (seconds)')
        title(['Events-based Seq #' num2str(iSeq)])
    end
    
end

all_seq.(type).K = K;
all_seq.(type).L = L;
all_seq.(type).lambda = lambda;
all_seq.(type).lambdaOrthoH = lambdaOrthoH;
all_seq.(type).lambdaOrthoW = lambdaOrthoW;
all_seq.(type).W = W_e;
all_seq.(type).H = H_e;
all_seq.(type).is_significant = is_significant;


%% Parts based version
fprintf('\n*************** Parts based ***************\n')
%parts based
rng(235); % fixed rng seed for reproduceability

lambdaOrthoH = 0;
lambdaOrthoW = 1; 
type = 'parts';

fprintf('Running seqNMF...K = %d  L = %d sec\n', K, L)

figure;
[W_p, H_p, ~,loadings,power]= seqNMF(X,'K',K,'L',Lneural,...
    'lambdaL1W', .1, 'lambda', lambda, 'maxiter', 100, 'showPlot', 1,...
    'lambdaOrthoH', lambdaOrthoH, 'lambdaOrthoW', lambdaOrthoW);

% test if any are significant
p = .05; % desired p value for factors

disp('Testing significance of factors on held-out data')
[pvals,is_significant] = test_significance(testNEURAL,W_p,p);

W_p = W_p(:,is_significant,:);
H_p = H_p(is_significant,:);
fprintf('Found %d/%d significant factors\n', sum(is_significant), length(is_significant))
% Plot the parts based results
if isempty(H_p)
    fprintf('<strong>No significant sequences detected.  Skipping...</strong>\n')
else
    figure(302);
    [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W_p(:,:,:),1);
    indSort = hybrid(:,3);
    tstart = 1; % plot data starting at this timebin
    WHPlot(W_p(indSort,:,:),H_p(:,tstart:end), X(indSort,tstart:end), ...
        0,trainPOS(:,floor(tstart*Fs/PosFs):end))
    title('Significant seqNMF factors, with raw data: Parts-based')
    
     saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_fac_parts' num2str(K) 'L' num2str(L) '_' type])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_fac_parts' num2str(K) 'L' num2str(L) '_' type '.png'])
    
        figure(303);
    WHPlot(W_p(indSort,:,:),H_p(:,tstart:end), ...
        helper.reconstruct(W_p(indSort,:,:),H_p(:,tstart:end)),...
        0,trainPOS(:,floor(tstart*Fs/PosFs):end))
    title('SeqNMF reconstruction: Parts-based')
    
         saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_recon_parts' num2str(K) 'L' num2str(L) '_' type])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_recon_parts' num2str(K) 'L' num2str(L) '_' type '.png'])
    
    figure(304)
     addpath(genpath([PARAMS.code_seqnmf_dir filesep 'misc_elm']));
    figure; HTriggeredSpec(H_p,trainPOS,Fs,Fs,ceil(L*Fs));
    title('Seq triggerred position')
    
    figure(305)
    figure; HTriggeredRaster(H_p,trainNEURAL(indSort,:),ceil(L*Fs));
    title('Seq triggered raster')
    
    % plot a reconstruction of the W matrix for sanity
    % plot the parts sequences with a bit more detail
    figure(306)
    for iSeq = 1:size(W_p,2)
        subplot(1, size(W_p,2),iSeq)
        % offset = 10^-1; % offset for each cell
        this_seq = squeeze(W_p(indSort,iSeq,:));
        imagesc((1:size(this_seq,1))/Fs,1:size(this_seq,2), this_seq);
        ylabel('sorted cell ID')
        xlabel('seq (seconds)');
        title(['Parts-based Seq #' num2str(iSeq)])
    end
    
    figure(307)
    tvec = (1:size(trainNEURAL,2))/Fs;
    temp_neural = trainNEURAL;
    temp_neural(temp_neural == 0) = NaN; 
    h = nan_imagesc_ec(temp_neural(indSort,:));
    x_tick = get(gca, 'xtick');
    set(gca, 'xticklabel', num2str(round(x_tick/Fs)'))
    hold on
    % get the H indices that exceed 2std
    [H_pks, H_pks_ind] =  sort(H_p(iSeq,:), 'descend');
    top_H = H_pks_ind(H_pks > (mean(H_p(iSeq,:)) +2*std(H_p(iSeq,:))));
    %    vline(tvec(top_H))
    %         plot(tvec(top_H), ones(1,length(top_H))*-10, '*', 'markersize', 20);
    r = rectangle('position', [tvec(top_H(1)), length(top_H)-10, L, 2],'edgecolor', PARAMS.red , 'facecolor', PARAMS.red);
    ylim([1-size(W_p,2)*2 size(trainNEURAL,1)])
end

% collect the output
all_seq.(type).K = K;
all_seq.(type).L = L;
all_seq.(type).lambda = lambda;
all_seq.(type).lambdaOrthoH = lambdaOrthoH;
all_seq.(type).lambdaOrthoW = lambdaOrthoW;
all_seq.(type).W = W_p;
all_seq.(type).H = H_p;
all_seq.(type).is_significant = is_significant;
all_seq.(type).indSort = indSort; 

save([PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'all_wake_seqs_' type], 'all_seq')

% clear this_pos this_data train* pos*
%% Get some REM data for this session


cd([REM_dir filesep subject filesep session])

% load the concatenated sleep data.
load('all_binary_post_REM.mat');
load('all_binary_post_SW.mat');
% restrict, transpose and rename the sleep types.
REM_data = all_binary_post_REM(:,keep_idx)';

if exist('idx', 'var')
    fprintf('Sorting data based on centroids in WAKE\n')
    REM_data = REM_data(idx,:);
end

SW_data = all_binary_post_SW(:,keep_idx)';

clear all_binary*

% plot the activity
figure(400)
imagesc(REM_data);
x_ticks = get(gca, 'XTick');
set(gca, 'Xticklabels', round(x_ticks/Fs))
xlabel('time (s)')
ylabel('cell id (unsorted)')
title('REM activity')
%% run Seq on REM and look for significant sequences using the same K and lambda from the task.

% this_REM = circshift(REM_data, floor(length(REM_data)/2),2);
this_REM = REM_data;
load("LFP_mats\all_t_post_REM.mat");
all_t_post_REM = smooth(all_t_post_REM, Fs)';

all_t_post_REM = [all_t_post_REM; all_t_post_REM]; % double for clear plotting.

splitN = floor(size(this_REM,2)*.75);

% make the test/trian sets.
trainNEURAL = this_REM(:,1:splitN);
testNEURAL = this_REM(:,(splitN+1):end);

trainLFP = all_t_post_REM(1:splitN);
testLFP = all_t_post_REM((splitN+1):end);

% Set some parameters
rng(235); % fixed rng seed for reproduceability
X = trainNEURAL;
%events vid
% lambdaOrthoH = .1; % favor events-based (these can take any value, don't need to be zero and one)
% lambdaOrthoW = 0;
% type = 'events';

%parts based
lambdaOrthoH = 0.1;
lambdaOrthoW = 0;
type = 'parts';

% Choose a K value.

% tic
% Ws = {};
% Hs = {};
% numfits = 3; %number of fits to compare
% for k = 1:4
%     display(sprintf('running seqNMF with K = %i',k))
%     for ii = 1:numfits
%         [Ws{ii,k},Hs{ii,k}] = seqNMF(X,'K',k, 'L', ceil(Ls(iL)*Fs),'lambda', 0,'maxiter',30,'showplot',0);
%         % note that max iter set low (30iter) for speed in demo (not recommended in practice)
%     end
%     inds = nchoosek(1:numfits,2);
%     for i = 1:size(inds,1) % consider using parfor for larger numfits
%             Diss(i,k) = helper.DISSX(Hs{inds(i,1),k},Ws{inds(i,1),k},Hs{inds(i,2),k},Ws{inds(i,2),k});
%     end
% end
% toc
% figure
% [~, best_k] = min(median(Diss))
% plot(1:4,Diss,'ko'), hold on
% h1 = plot(1:4,median(Diss,1),'k-','linewidth',2);
% h2 = plot([best_k,best_k],[0,0.5],'r--');
% legend([h1 h2], {'median Diss','true K'})
% xlabel('K')
% ylabel('Diss')
%%
% Run seqNMF using the best K +1

Ls = 2:2:16;
if ~exist("best_k", 'var')
    best_k = 2;
end

for iL = 1:length(Ls)
    fprintf('Running seqNMF...K = %d  L = %d sec\n', best_k+1, Ls(iL))
    [W, H, ~,loadings,~]= seqNMF(X,'K',best_k+1,'L',ceil(Ls(iL)*Fs),...
        'lambdaL1W', .1, 'lambda', lambda, 'maxiter', 100, 'showPlot', 0,...
        'lambdaOrthoH', lambdaOrthoH, 'lambdaOrthoW', lambdaOrthoW);
    
    % test if any are significant
    p = .05; % desired p value for factors
    
    disp('Testing significance of factors on held-out data')
    [pvals,is_significant] = test_significance(testNEURAL,W,p);
    
    W = W(:,is_significant,:);
    H = H(is_significant,:);
    fprintf('Found %d/%d significant factors\n', sum(is_significant), length(is_significant))
    all_sweeps_out{iL}.K = K;
    all_sweeps_out{iL}.L = Ls(iL);
    all_sweeps_out{iL}.W = W;
    all_sweeps_out{iL}.H = H;
    all_sweeps_out{iL}.indSort = indSort;
    all_sweeps_out{iL}.loadings = loadings;
    all_sweeps_out{iL}.sig = is_significant;
    all_sweeps_out{iL}.Train = trainNEURAL;
    all_sweeps_out{iL}.Trest = testNEURAL;
    
end


all_sweeps.(type) = all_sweeps_out; % workaround for parfor

sig_Ls = []; 
for iSeq = 1:length(all_sweeps.(type))
    if sum(all_sweeps.(type){iSeq}.sig) >0
        sig_Ls = [sig_Ls iSeq]; 
        fprintf('Significant Seq (%d) found using k: %d L: %d\n', sum(all_sweeps.(type){iSeq}.sig), all_sweeps.(type){iSeq}.K, all_sweeps.(type){iSeq}.L)
    end
end
%%
% Reconstruct a particular L sweep.

for iL = sig_Ls % pick an L to work with
    
    fprintf('<strong>Plotting %d/%d significant factors for K = %d L = %ds</strong>\n', sum(all_sweeps.(type){iL}.sig), length(all_sweeps.(type){iL}.sig),all_sweeps.(type){iL}.K, all_sweeps.(type){iL}.L )
    
    if sum(all_sweeps.(type){iL}.sig) >0
        
        figure;
        [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(all_sweeps.(type){iL}.W(:,:,:),1);
        indSort = hybrid(:,3);
        tstart = 1; % plot data starting at this timebin
        WHPlot(all_sweeps.(type){iL}.W(indSort,:,:),all_sweeps.(type){iL}.H(:,tstart:end), all_sweeps.(type){iL}.Train(indSort,tstart:end), ...
            0)
        title(['REM Significant seqNMF factors, with raw data: ' type '-based L:' num2str(Ls(iL)) 's'])
        mkdir([PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF'])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'REM_Seq_factor_parts_K' num2str(K) 'L' num2str(Ls(iL))])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'REM_Seq_factor_parts_K' num2str(K) 'L' num2str(Ls(iL)) '.png'])
        
        figure;
        WHPlot(all_sweeps.(type){iL}.W(indSort,:,:),all_sweeps.(type){iL}.H(:,tstart:end), ...
            helper.reconstruct(all_sweeps.(type){iL}.W(indSort,:,:),all_sweeps.(type){iL}.H(:,tstart:end)),...
            0)
        title(['REM SeqNMF reconstruction: ' type '-based L:' num2str(Ls(iL)) 's'])
        
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'REM_Seq_recon_parts_K' num2str(K) 'L' num2str(Ls(iL))])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'REM_Seq_recon_parts_K' num2str(K) 'L' num2str(Ls(iL)) '.png'])
        
        % close all
    end
end
% reconstruct the data using color coded anx and safety cells.

for iSeq = 1:size(all_sweeps.(type){iL}.W(indSort,:,:),2)
    figure(1000+iSeq)
    subplot(5,4,1)
    text(.1,.8, {[type '-based']});
     text(.1,.4, ['L: ' num2str(Ls(iL)) 's'])
          text(.1,0, ['REM Seq # ' num2str(iSeq)])
    axis off; 
    
    ax1(1) =  subplot(5,4,2:4);
    imagesc((1:size(REM_data, 2))/Fs,1,all_t_post_REM(1,:))
    set(gca, 'xtick', [], 'ytick', []); 
    title('Theta power')
    colormap('jet')

    ax2(1) = subplot(5,4,[5 9 13 17]);
    this_seq = squeeze(all_sweeps.(type){iL}.W(indSort,iSeq,:));
    this_seq(this_seq > 0.1) = 1; 
    imagesc((1:size(this_seq, 2))/Fs,1:size(this_seq, 1), this_seq.*([1:size(REM_data,1)]+(floor(size(this_seq,1)/2)))');
    ylabel('cell ID (sorted)')
    xlabel('seq time (s)');
    
   ax1(2) = subplot(5,4,[6:8 10:12 14:16 18:20]);
    imagesc((1:size(REM_data, 2))/Fs,1:size(REM_data, 1), REM_data(indSort,:).*([1:size(REM_data,1)]+(floor(size(this_seq,1)/2)))');
    set(gca, 'ytick', []); 
    xlabel('concatenated time (s)')
    cfg_plot.ft_size = 12; 
    SetFigure(cfg_plot, gcf);
    set(gcf, 'position', [680 500 1127 480])
    
        linkaxes(ax1, 'x')
    
    colormap('jet');
    
    colormap(ax2(1), 'parula');
    colormap(ax1(2), 'parula');
    
%        linkaxes(ax2, 'x')
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'REM_Seq_raw_parts_K' num2str(K) 'L' num2str(Ls(iL)) '_' type])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'REM_Seq_raw_parts_K' num2str(K) 'L' num2str(Ls(iL)) '_' type '.png'])
end


%% same plot but ordered based on sequences in wake

sort_REM = REM_data(all_seq.parts.indSort,:); 

    figure(2000)
    subplot(5,4,1)
    text(.1,.8, {[type '-based']});
     text(.1,.4, ['L: ' num2str(Ls(iL)) 's'])
          text(.1,0, ['REM Seq # ' num2str(iSeq) '_Wake_Sort'])
    axis off; 
    
    ax1(1) =  subplot(5,4,2:4);
    imagesc((1:size(sort_REM, 2))/Fs,1,all_t_post_REM(1,:))
    set(gca, 'xtick', [], 'ytick', []); 
    title('Theta power')
    colormap('jet')

    ax2(1) = subplot(5,4,[5 9 13 17]);
    this_seq = squeeze(all_sweeps.(type){iL}.W(all_seq.parts.indSort,iSeq,:));
    this_seq(this_seq > 0.1) = 1; 
    imagesc((1:size(this_seq, 2))/Fs,1:size(this_seq, 1), this_seq.*([1:size(REM_data,1)]+(floor(size(this_seq,1)/2)))');
    ylabel('cell ID (sorted)')
    xlabel('seq time (s)');
    
   ax1(2) = subplot(5,4,[6:8 10:12 14:16 18:20]);
    imagesc((1:size(REM_data, 2))/Fs,1:size(REM_data, 1), REM_data(indSort,:).*([1:size(REM_data,1)]+(floor(size(this_seq,1)/2)))');
    set(gca, 'ytick', []); 
    xlabel('concatenated time (s)')
    cfg_plot.ft_size = 12; 
    SetFigure(cfg_plot, gcf);
    set(gcf, 'position', [680 500 1127 480])
    
        linkaxes(ax1, 'x')
    
    colormap('jet');
    
    colormap(ax2(1), 'parula');
    colormap(ax1(2), 'parula');
    
%        linkaxes(ax2, 'x')
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'REM_wake_sort_Seq_raw_parts_K' num2str(K) 'L' num2str(Ls(iL)) '_' type])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'REM_wake_sort_Seq_raw_parts_K' num2str(K) 'L' num2str(Ls(iL)) '_' type '.png'])


% factor plots Not sure if useful for REM and LFP power.
% % Plot factor-triggered song examples and rastors
% figure(503)
% addpath(genpath([PARAMS.code_seqnmf_dir filesep 'misc_elm']));
% figure; HTriggeredSpec(all_sweeps_parts{iL}.H,trainLFP,Fs,Fs,ceil(all_sweeps_parts.L*Fs)); % spectrogram (just LFP power in this case)
% title('Factor Triggered LFP power')
% figure;
% HTriggeredRaster(all_sweeps_parts{iL}.H,trainNEURAL(indSort,:),ceil(all_sweeps_parts.L*Fs)); % raster
% title('Factor Triggered raster')
% clearvars -except all_sweeps trainNEURAL testNEURAL Ls
%%

if cfg.score
% Check the Sequence Score
tic
nRepsShuff = cfg.score_nShuf;
nRepsColShuff = cfg.score_nShuf;  % 15 for just making an estimate.  100 for test?

nIter = 50;

L = all_sweeps.(type){iL}.L; % from earlier run.
K = all_sweeps.(type){iL}.K;

X = this_REM;
[N T] = size(X);

% do seqNMF
tmp = [];

parfor iteri = 1:nIter
    rng('shuffle')
    [~, ~, ~,~,tmp(iteri)] = seqNMF(X, 'L', L, 'K', K, 'lambda', 0, 'showPlot',0);
    fprintf('part 1 inter # %d\n', iteri)
end

PEx = max(tmp);


% do seqNMF on shuffled data
PExShuff = [];
parfor repi = 1:nRepsShuff
    Xshuff = [];
    for ni = 1:N
        timeshuff = randperm(T);
        Xshuff(ni,:) = X(ni, timeshuff);
    end
    tmp = [];
    
    for iteri = 1:nIter
        rng('shuffle')
        [~, ~, ~,~,tmp(iteri)] = seqNMF(Xshuff, 'L', L, 'K', K, 'lambda', 0.0, 'showPlot',0);
        fprintf('part 2 inter # %d\n', iteri)
    end
    PExShuff(repi) = max(tmp);
    fprintf('part 2 rep # %d\n', repi)
    
end

% do seqNMF on col shuffled data
PExColShuff = [];
parfor repi = 1:nRepsColShuff
    Xshuff = X(:,[1:L (L + randperm(T-L))]); % don't shuffle to first L bins... these cannot be explained by seqNMF
    tmp = [];
    
    for iteri = 1:nIter
        rng('shuffle')
        [~, ~, ~,~,tmp(iteri)] = seqNMF(Xshuff, 'L', L, 'K', K, 'lambda', 0.0, 'showPlot',0);
        fprintf('part 3 inter # %d\n', iteri)
        
    end
    PExColShuff(repi) = max(tmp);
    fprintf('part 3 rep # %d\n', repi)
    
end
NoiseFloor = median(PExShuff);
SyncFloor = median(PExColShuff);
PAS = (PEx-SyncFloor)./...
    (PEx-NoiseFloor)
fprintf('Score: %2.4f\n', PAS)
toc
all_sweeps.(type){iL}.PAS = PAS;
end
%% save the output

save([PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'all_REM_sweeps_' type ], 'all_sweeps')
