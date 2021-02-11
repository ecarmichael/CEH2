%% SeqNMFE sandbox: Jisu Sleep Seq Analyses

close all
restoredefaultpath
global PARAMS  % these are global parameters that can be called into any function.  I limit these to directories for storing, loading, and saving files and codebases.
os = computer;

if ismac
    PARAMS.data_dir = '/Users/jericcarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/Users/jericcarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    
elseif strcmp(os, 'GLNXA64')
    
    %     PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/7_12_2019_PV1069_LTD5'; % where to find the raw data
    %     PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/EV/';
    PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/'; % raw data location.
    PARAMS.inter_dir = '/home/ecarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/home/ecarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    
else
    PARAMS.data_dir = 'D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\Processed place'; % where to find the raw data
    PARAMS.raw_data_dir = 'D:\Dropbox (Williams Lab)\Jisoo\JisooProject2020\RawData\pv1060'; % raw data location.
    PARAMS.inter_dir = 'D:\Dropbox (Williams Lab)\Jisoo\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter'; % where to put intermediate files
    PARAMS.stats_dir = 'D:\Dropbox (Williams Lab)\Jisoo\Upload\pv1002\Inter\Stats'; % where to put the statistical output .txt
    PARAMS.code_base_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2'; % where the multisite repo can be found
    PARAMS.code_seqnmf_dir = 'C:\Users\ecarm\Documents\GitHub\seqNMF'; % where the multisite repo can be found
    
end

% colours
PARAMS.L_grey = [0.8 0.8 0.8];
PARAMS.D_grey = [0.2 0.2 0.2];
PARAMS.blue = [0.3639    0.5755    0.7484];
PARAMS.red = [0.9153    0.2816    0.2878];
PARAMS.green= [0.4416    0.7490    0.4322];
PARAMS.gold = [1.0000    0.5984    0.2000];

rng(11,'twister') % for reproducibility

% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_CEH2_dir));
cd(PARAMS.raw_data_dir) % move to the data folder

clear d os
%% load a nice session

%%%%%% type of analysis to run %%%%%%%
% for place cells only
this_analysis = 'place';

% for anxiety only
this_analysis = 'anxiety';% default. comment out if you want place only. 


% this_sess = 'D:\Dropbox (Williams Lab)\Jisoo\Upload\pv1002\3_26_2019';
% this_sess = 'D:\Dropbox (Williams Lab)\Jisoo\Upload\pv1043\6_11_2019_PV1043_LTD1\LTD1'
% this_sess = 'D:\Dropbox (Williams Lab)\Jisoo\JisooProject2020\RawData\pv1043\6_15_2019_PV1043_LTD5\H13_M30_S2_LTD5';
% this_sess = 'D:\Dropbox (Williams Lab)\Jisoo\JisooProject2020\RawData\pv1043\6_13_2019_PV1043_LTD3\H13_M18_S10_LTD3';
% this_sess = 'D:\Dropbox (Williams Lab)\Jisoo\JisooProject2020\RawData\pv1060\11_23_2019_PV1060_HATD5\H13_M25_S1_HATD5';
% this_sess = 'D:\Dropbox (Williams Lab)\Jisoo\JisooProject2020\RawData\pv1069\10_22_2019_PV1069_HATSwitch\H13_M4_S44_HATD6';
this_sess = 'D:\Dropbox (Williams Lab)\Jisoo\JisooProject2020\RawData\pv1060\11_26_2019_PV1060_HATSwitch\H13_M5_S15_HATSwitch';

dir_parts = strsplit(this_sess, filesep);
task = strsplit(dir_parts{end}, '_');
task = task{end}; 
if contains(task, 'HATD6')
    task = 'HATDSwitch';
end
subject = dir_parts{end-2};
session = dir_parts{end-1};

% processed data dir
% cd('D:\Dropbox (Williams Lab)\Jisoo\JisooProject2020\2020_Results_aftercutting\4.PlaceCell\pv1043\LTD5');
% cd('D:\Dropbox (Williams Lab)\Jisoo\JisooProject2020\2020_Results_aftercutting\4.PlaceCell\pv1043\LTD3');
this_process_dir ='D:\Dropbox (Williams Lab)\Jisoo\JisooProject2020\2020_Results_aftercutting\';

% REM file dir
this_rem_dir = ['D:\Dropbox (Williams Lab)\Jisoo\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter' filesep upper(subject) filesep session];


% load the main data
cd(this_sess);

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

cd(this_process_dir)

if contains(this_sess, 'HAT')
    %     if contains(this_sess, 'Switch')
    %         session  = strrep(session, 'HATSwitch', 'HATDSwitch'); % switch naming convention if needed for 'Switch'
    %     end
    %     sess_name = strsplit(session, '_');
    %     S = dir(fullfile(this_process_dir,[lower(sess_name{end-1}) '_' sess_name{end} '*cell.mat']));
    %     load(fullfile(this_process_dir,S.name));
    %         fprintf('Using a HAT session. Loading cell file: <strong>%s</strong>\n', S.name);
    if exist([this_process_dir  '3.Anxiety_SafetyCell' filesep subject filesep 'Anxiety_output' filesep 'STD1'], 'dir')
        load([this_process_dir '3.Anxiety_SafetyCell' filesep subject filesep 'Anxiety_output' filesep 'STD2' filesep 'Anxiety_Output_std2.mat'])
    else
        load([this_process_dir  '3.Anxiety_SafetyCell' filesep subject filesep 'Anxiety_output'  filesep 'Anxiety_Output_std2.mat'])
    end
    % get place cells as well.
    
    if contains(task, 'HATD')% workaround for naming of HAT + D
        load([this_process_dir filesep '4.PlaceCell' filesep subject filesep task filesep 'spatial_analysis.mat'])
        load([this_process_dir filesep '4.PlaceCell' filesep subject filesep task filesep 'SA.mat'])
    else
        load([this_process_dir filesep '4.PlaceCell' filesep subject filesep strrep(task, 'HAT', 'HATD') filesep 'spatial_analysis.mat'])
        load([this_process_dir filesep '4.PlaceCell' filesep subject filesep strrep(task, 'HAT', 'HATD') filesep 'SA.mat'])
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
if contains(this_sess, 'HAT'); keep_idx = (Anxiety_Output.AnxietyCell_day_switch); end

if contains(this_analysis, 'place')
    keep_idx = place_cell_idx; 
end
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

figure(111)
hold on
% plot(ms.time, pos, 'color', [.2 .2 .2])
plot(L_laps_time/1000, nan(length(L_laps_pos),1), 'b')
plot(R_laps_time/1000, nan(length(R_laps_pos),1), 'r')

plot(L_laps_time/1000, L_laps_pos, '.b')
plot(R_laps_time/1000, R_laps_pos, '.r')
xlabel('time (s)')
ylabel('position')
legend({'Left Laps', 'Right Laps'})
title('Position (L/R laps)')

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
figure(112);
plot_pos_mat = pos_mat;
plot_pos_mat(plot_pos_mat == 0) =NaN;
% plot_pos_mat(~isnan(plot_pos_mat)) = 1;
nan_imagesc_ec(plot_pos_mat);
axis xy
x_ticks = get(gca, 'xtick');
set(gca, 'xticklabels', x_ticks/Fs)
set(gca, 'color', 'w')
colormap('jet')
title('position matrix')
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

all_seq.subject = subject;
all_seq.session = session;
all_seq.task = dir_parts{end}(end-3:end);
all_seq.analysis_date = datetime;
all_seq.data = this_data;
all_seq.pos = this_pos;

%% Procedure for choosing lambda

% nLambdas = 30; % increase if you're patient
% K = 3;
% X = trainNEURAL;
% lambdas = sort([logspace(-2,-4,nLambdas)], 'ascend');
% loadings = [];
% regularization = [];
% cost = [];
%     [N,T] = size(X);
%
% parfor li = 1:length(lambdas)
%     [W, H, ~,loadings(li,:),~]= seqNMF(X,'K',K,'L',Lneural,...
%         'lambdaL1W', .1, 'lambda', lambdas(li), 'maxiter', 100, 'showPlot', 0);
%     [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X,W,H);
%     display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
% end
%% plot Lambda cost

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
%
% % get the lambda point after the cross point.
% lambda_diff = Rs - Cs;
% idx = find(lambda_diff < 0);
%
% best_lambda = lambdas(idx(1));


%% Events based version

lambda = 0.001;% from the above.  "Go just above the cross point"
fprintf('\n*************** Events based ***************\n')
type = 'events';
% event based
lambdaOrthoH = .1; % favor events-based (these can take any value, don't need to be zero and one)
lambdaOrthoW = 0;
fprintf('Running seqNMF...K = %d  L = %d sec\n', K, L)
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
    
    figure;
    WHPlot(W_e(indSort,:,:),H_e(:,tstart:end), ...
        helper.reconstruct(W_e(indSort,:,:),H_e(:,tstart:end)),...
        0,trainPOS(:,floor(tstart*Fs/PosFs):end))
    title('SeqNMF reconstruction: Event-based')
    
    
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
    
    
        figure(303);
    WHPlot(W_p(indSort,:,:),H_p(:,tstart:end), ...
        helper.reconstruct(W_p(indSort,:,:),H_p(:,tstart:end)),...
        0,trainPOS(:,floor(tstart*Fs/PosFs):end))
    title('SeqNMF reconstruction: Parts-based')
    
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
    h.XData = [tvec, tvec(end)+Fs]; 
%     h = imagesc(tvec, 1:size(trainNEURAL,1), (trainNEURAL(indSort,:)));
        set(gca, 'color', 'w')

% x_ticks = get(gca, 'xtick');
% set(gca, 'xticklabels', num2str(floor(x_ticks/Fs)))
    hold on
    % get the H indices that exceed 2std
    [H_pks, H_pks_ind] =  sort(H_p(iSeq,:), 'descend');
    top_H = H_pks_ind(H_pks > (mean(H_p(iSeq,:)) +2*std(H_p(iSeq,:))));
    %    vline(tvec(top_H))
    %         plot(tvec(top_H), ones(1,length(top_H))*-10, '*', 'markersize', 20);
    r = rectangle('position', [tvec(top_H(1)), -10/iSeq, L, 10/size(H_p,3)],'edgecolor', PARAMS.red , 'facecolor', PARAMS.red);
    ylim([-10 size(trainNEURAL,1)])
    colormap([1 1 1; 0 0 0])
    set(gcf, 'position', [303   562   964   350]); 
%     colormap(flipud(gray))


%% re-order the raw data. 

for iSeq = 1:size(W_p(indSort,:,:),2)
    figure(310+iSeq)
    subplot(5,4,1)
    text(.1,.8, {[type '-based']}, 'fontsize', 10, 'fontweight', 'bold');
     text(.1,.4, ['L: ' num2str(L) 's'], 'fontsize', 10, 'fontweight', 'bold')
          text(.1,0, ['Seq # ' num2str(iSeq)], 'fontsize', 10, 'fontweight', 'bold')
    axis off; 
    
     ax(2) = subplot(5,4,2:4);
     temp_pos = this_pos;
     temp_pos(temp_pos == 0) = NaN; 
    h = nan_imagesc_ec(temp_pos);
    h.XData = (1:length(this_pos)+1)/Fs;
    axis xy
    x_ticks = get(gca, 'xtick');
    set(gca, 'xticklabels', x_ticks/Fs)
    set(gca, 'color', 'w')
    colormap('jet')
%     imagesc((1:size(this_pos, 2))/Fs,1:size(this_pos,1),this_pos.*([1:size(this_pos,2)]+500))
    set(gca, 'xtick', [], 'ytick', []); 
    title('position')

    subplot(5,4,[5 9 13 17])
    this_seq = squeeze(W_p(indSort,iSeq,:));
    this_seq(this_seq > 0.1) = 1; 
    imagesc((1:size(this_seq, 2))/Fs,1:size(this_seq, 1), this_seq.*([1:size(this_seq,1)]+(floor(size(this_seq,1)/2)))');
    ylabel('cell ID (sorted)')
    xlabel('seq time (s)');
    
    ax(3) = subplot(5,4,[6:8 10:12 14:16 18:20]);
    this_seq = this_data(indSort,:);
    imagesc((1:size(this_seq, 2))/Fs,1:size(this_seq, 1), this_seq.*([1:size(this_seq,1)]+(floor(size(this_seq,1)/2)))');
    set(gca, 'ytick', []); 
    xlabel('concatenated time (s)')
    cfg_plot.ft_size = 12; 
    SetFigure(cfg_plot, gcf);
    set(gcf, 'position', [680 500 1127 480])
    
    linkaxes(ax, 'x')
end

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

mkdir([PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF'])
save([PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'all_wake_seqs_' this_analysis], 'all_seq')

% clear this_pos this_data train* pos*
%% Get some REM data for this session


cd(this_rem_dir)

% load the concatenated sleep data.
load('all_binary_post_REM.mat');
load('all_binary_post_SW.mat');
% restrict, transpose and rename the sleep types.
REM_data = all_binary_post_REM(:,keep_idx)';

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

Ls = 6:2:16;
if ~exist("best_k", 'var')
    best_k = 1;
end

parfor iL = 1:length(Ls)
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
    all_sweeps_out{iL}.loadings = loadings;
    all_sweeps_out{iL}.sig = is_significant;
    all_sweeps_out{iL}.Train = trainNEURAL;
    all_sweeps_out{iL}.Trest = testNEURAL;
    all_sweeps_out{iL}.lambdaOrthoH = lambdaOrthoH;
    all_sweeps_out{iL}.lambdaOrthoW = lambdaOrthoW;

end

all_sweeps.(type) = all_sweeps_out; % workaround for parfor

for iSeq = 1:length(all_sweeps.(type))
    if sum(all_sweeps.(type){iSeq}.sig) >0
        fprintf('Significant Seq (%d) found using k: %d L: %d\n', sum(all_sweeps.(type){iSeq}.sig), all_sweeps.(type){iSeq}.K, all_sweeps.(type){iSeq}.L)
    end
end
%%
% Reconstruct a particular L sweep.

for iL = 3 % pick an L to work with
    
    fprintf('<strong>Plotting %d/%d significant factors for K = %d L = %ds</strong>\n', sum(all_sweeps.(type){iL}.sig), length(all_sweeps.(type){iL}.sig),all_sweeps.(type){iL}.K, all_sweeps.(type){iL}.L )
    
    if sum(all_sweeps.(type){iL}.sig) >0
        
        figure(801);
        [max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(all_sweeps.(type){iL}.W(:,:,:),1);
        indSort = hybrid(:,3);
        tstart = 1; % plot data starting at this timebin
        WHPlot(all_sweeps.(type){iL}.W(indSort,:,:),all_sweeps.(type){iL}.H(:,tstart:end), all_sweeps.(type){iL}.Train(indSort,tstart:end), ...
            0)
        title(['REM Significant seqNMF factors, with raw data: ' type '-based L:' num2str(Ls(iL)) 's'])
        mkdir([PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF'])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_factor_parts_K' num2str(K) 'L' num2str(Ls(iL)) '_' this_analysis])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_factor_parts_K' num2str(K) 'L' num2str(Ls(iL)) '_' this_analysis '.png'])
        
        figure(802);
        WHPlot(all_sweeps.(type){iL}.W(indSort,:,:),all_sweeps.(type){iL}.H(:,tstart:end), ...
            helper.reconstruct(all_sweeps.(type){iL}.W(indSort,:,:),all_sweeps.(type){iL}.H(:,tstart:end)),...
            0)
        title(['REM SeqNMF reconstruction: ' type '-based L:' num2str(Ls(iL)) 's'])
        
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_recon_parts_K' num2str(K) 'L' num2str(Ls(iL)) '_' this_analysis])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_recon_parts_K' num2str(K) 'L' num2str(Ls(iL)) '_' this_analysis '.png'])
        
        % close all
    end
end
% reconstruct the data using color coded anx and safety cells.

for iSeq = 1:size(all_sweeps.(type){iL}.W(indSort,:,:),2)
    figure(1000+iSeq)
    subplot(5,4,1)
    text(.1,.8, {[type '-based']});
     text(.1,.4, ['L: ' num2str(iL) 's'])
          text(.1,0, ['Seq # ' num2str(iSeq)])
    axis off; 
    
    ax2(1) =  subplot(5,4,2:4);
    imagesc((1:size(REM_data, 2))/Fs,1,all_t_post_REM(1,:))
    set(gca, 'xtick', [], 'ytick', []); 
    title('Theta power')
    colormap('jet')

    subplot(5,4,[5 9 13 17])
    this_seq = squeeze(all_sweeps.(type){iL}.W(indSort,iSeq,:));
    this_seq(this_seq > 0.1) = 1; 
    imagesc((1:size(this_seq, 2))/Fs,1:size(this_seq, 1), this_seq.*([1:size(REM_data,1)]+(floor(size(this_seq,1)/2)))');
    ylabel('cell ID (sorted)')
    xlabel('seq time (s)');
    
   ax2(2) = subplot(5,4,[6:8 10:12 14:16 18:20]);
    imagesc((1:size(REM_data, 2))/Fs,1:size(REM_data, 1), REM_data(indSort,:).*([1:size(REM_data,1)]+(floor(size(this_seq,1)/2)))');
    set(gca, 'ytick', []); 
    xlabel('concatenated time (s)')
    cfg_plot.ft_size = 12; 
    SetFigure(cfg_plot, gcf);
    set(gcf, 'position', [680 500 1127 480])
    
%        linkaxes(ax2, 'x')
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_raw_parts_K' num2str(K) 'L' num2str(Ls(iL)) '_' this_analysis])
        saveas(gcf, [PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'Seq_raw_parts_K' num2str(K) 'L' num2str(Ls(iL)) '_' this_analysis '.png'])
        
end

%     figure(809)
%     hold on
%     c_ord = parula(size(REM_data,1));
%     REM_tvec = (1:size(REM_data, 2))/Fs;
%     REM_sort = REM_data(indSort,:);
%     REM_sort(REM_sort == 0) = NaN;
% for ii = size(REM_data,1):-1:1
%     plot(REM_tvec, (REM_sort(ii, :)*.8)-ii, 'color', c_ord(ii,:), 'linewidth', 10)
% end
% xlim([min(REM_tvec) max(REM_tvec)]);
% ylim([ -size(REM_data,1) 0 ]);
% axis xy

ylabel('sorted cell ID')
xlabel('seq frames');
title(type)

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
% Check the Sequence Score
tic
nRepsShuff = 15;
nRepsColShuff = 15;  % 15 for just making an estimate.  100 for test?

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
    (PEx-NoiseFloor);
fprintf('Score: %2.4f\n', PAS)
toc
all_sweeps.(type){iL}.PAS = PAS;
%% save the output

save([PARAMS.inter_dir filesep subject filesep session filesep 'SeqNMF' filesep 'all_REM_sweeps_' this_analysis], 'all_seq')