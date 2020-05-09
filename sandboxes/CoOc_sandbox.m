%% MS_LFP_event_reactivation_sandbox

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
    PARAMS.raw_data_dir = 'D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\LT&sleep'; % raw data location.
    PARAMS.csc_data_dir = 'D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\LFP'; % where are the LFP files. If blank will look in the same folder as raw_data.
    PARAMS.inter_dir = 'D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\Inter'; % where to put intermediate files
    PARAMS.stats_dir = 'D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\Inter\Stats'; % where to put the statistical output .txt
    PARAMS.code_base_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2'; % where the multisite repo can be found
    PARAMS.code_seqnmf_dir = 'C:\Users\ecarm\Documents\GitHub\seqNMF'; % where the multisite repo can be found
    
end


rng(11,'twister') % for reproducibility


% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_CEH2_dir));
cd(PARAMS.raw_data_dir) % move to the data folder

% try the newer NLX loaders for UNIX
[~, d] = version;
if str2double(d(end-3:end)) >2014 && strcmp(os, 'GLNXA64')
    rmpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared/io/neuralynx')
    addpath(genpath('/Users/jericcarmichael/Documents/NLX_loaders_UNIX_2015'))
    disp('Version is greater than 2014b on UNIX so use updated loaders found here:')
    which Nlx2MatCSC
end

clear d os

%% load a nice session

load('D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\Inter\537\12_5_2019_537day1\ms_resize.mat');

% load('D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\Inter\537\12_5_2019_537day1\ms_trk.mat');
%%
for iB = 1:length(ms_seg_resize.RawTraces)
    len = ms_seg_resize.NLX_csc{iB}.tvec(end) - ms_seg_resize.NLX_csc{iB}.tvec(1);
    str_len = length(strcat(ms_seg_resize.file_names{iB},ms_seg_resize.pre_post{iB}, ms_seg_resize.hypno_label{iB}));
    fprintf('<strong>%s %s - %s:</strong>', ms_seg_resize.file_names{iB},ms_seg_resize.pre_post{iB}, ms_seg_resize.hypno_label{iB})
    fprintf(repmat(' ', 1,abs(str_len-23)))
    fprintf('SWD:%3d  %2.1f/s   SWR:%3d  %2.1f/s   low_G:%3d  %2.1f/s \n',...
        length(ms_seg_resize.SWD_evts{iB}.tstart),length(ms_seg_resize.SWD_evts{iB}.tstart)/len,....
        length(ms_seg_resize.SWR_evts{iB}.tstart),length(ms_seg_resize.SWR_evts{iB}.tstart)/len,...
        length(ms_seg_resize.low_gamma_evts{iB}.tstart),length(ms_seg_resize.low_gamma_evts{iB}.tstart)/len)
end
fprintf('_______________________________________________________________________________\n')
%% start with one block.
iSeg = 1;

keep_idx = 1:size(ms_seg_resize.RawTraces,1); % actually this is a remove index
keep_idx =keep_idx(find((keep_idx ~= iSeg)));

cfg_rem = [];
ms_seg = MS_remove_data_sandbox(cfg_rem, ms_seg_resize, keep_idx);

ms_seg = MS_de_cell(ms_seg);

% binarize the trace

ms_seg = msExtractBinary_detrendTraces(ms_seg);

% check for inactive cells and remove from ms.SFPs just using sum of
% binary > 0;

cfg_SFP = [];
cfg_SFP.fnc = '==';
cfg_SFP.remove_val = 0;
ms_seg = MS_update_SFP(cfg_SFP, ms_seg);

csc = ms_seg.NLX_csc;
SWRs = ms_seg.SWR_evts;
SWRs.tstart = SWRs.tstart - csc.tvec(1);
SWRs.tend = SWRs.tend - csc.tvec(1);
nlx_evts = ms_seg.NLX_evt;
nlx_evts.t{end} = nlx_evts.t{end} - csc.tvec(1);
% correc the tvec;
csc.tvec = csc.tvec - csc.tvec(1);


Ca_TS = MS_Binary2TS([], ms_seg);
Ca_TS.usr = [];
% make a subset of neurons (for speed)
nS = 951;
Ca_1 = Ca_TS; Ca_1.t = []; Ca_1.label = [];
for iC  = nS:-1:1
    Ca_1.t{iC} = Ca_TS.t{iC} -csc.tvec(1) ; Ca_1.label{iC} = Ca_TS.label{iC};
end



%% try some stuff
figure(12)
% little checker.
traces = 1:nS;
ax(1)= subplot(4,1,1);
cfg.target = 'LFP';
PlotTSDfromIV(cfg, SWRs, csc);
%     plot(csc.tvec - csc.tvec(1), csc.data(2,:))
xlim([csc.tvec(1) csc.tvec(end)])


ax(3)= subplot(4,1,2:4);

for iT  = traces
    hold on
    this_B = ms_seg.Binary(:, iT);
    this_B(this_B == 0) = NaN; 
    plot(ms_seg.time/1000,    (this_B*(0.01*length(traces)))+iT, 'linewidth', 2)
    this_B =[]; 
end
xlim([ms_seg.time(1)/1000, ms_seg.time(end)/1000]);
ylim([traces(1) traces(end)]);
linkaxes(ax, 'x')

%% event by event plots w/checks

% [SWRs, removed_idx] = MS_check_IV([], csc, ms_seg, SWRs);

remove_idx = [22    30    67    68    69    90   132   133   141   142]; % save time for 537-12-5-2019 day 1
SWRs.tstart(remove_idx) = [];
SWRs.tend(remove_idx) = [];
%% get some basic stats for active cells per SWR;
cfg = [];
cfg.t_win = [-0.1 0.5]; 
cfg.plot = 1;
spike_counts = MS_event_hist(cfg, Ca_TS, csc, SWRs);

t_edges = cfg.t_win(1):0.01:cfg.t_win(end);
t_centers = t_edges(1:end-1)+0.01/2;
bar(t_centers*1000,  mean(sum(spike_counts,3))*100);
vline(0)
xlabel('time (ms)')
ylabel('mean % active cells')

%% get shuffle distribution {remember to turn set cfg.plot = 0. }. Can parfor. parfor can do 100 shuffles in 12s with 8 workers and 1000 shuffles in 124s. regular for loop takes 85s. 
cfg.plot = 0; 
nShuffle = 100; 
shuf_counts = nan(nShuffle, length(t_centers));
tic
tvec = csc.tvec; 
parfor iShuffle = 1:nShuffle
    shuf_IV = MS_get_random_epochs(tvec, SWRs);
    shuf_counts(iShuffle,:) = mean(sum(MS_event_hist(cfg, Ca_TS, csc, shuf_IV),3));
    
end
toc

hold on
bar(t_centers*1000, mean(shuf_counts)*100);
legend('SWR', '1000 shuffle'); 
%%
swr_centers = IVcenters(SWRs); % convert to centered events;


% make a wider SWRs IV

SWRs_wide = SWRs;
SWRs_wide.tend = SWRs.tend + 1;

% get the time idx that matches the SWR centers (use this if you just want
% one frame before and one after or something.
swr_ms_idx_centers = nearest_idx3(swr_centers, nlx_evts.t{end});

% alternative:
% get the idx for the start and end of the event.
swr_ms_idx_tstart = nearest_idx3(SWRs.tstart, nlx_evts.t{end});

swr_ms_idx_tend = nearest_idx3(SWRs.tend, nlx_evts.t{end});

% initialize some matricies to store the co-activity.
co_mat = NaN(size(ms_seg.Binary,2),size(ms_seg.Binary,2),length(swr_ms_idx_tstart)); % Make an empty matrix for co-activity
corr_mat = NaN(size(ms_seg.Binary,2),size(ms_seg.Binary,2),length(swr_ms_idx_tstart)); % Make an empty matrix for correlations
Q_mat = NaN(size(ms_seg.Binary,2),length(swr_ms_idx_tstart));
% %% quick check for Q mat
% 
% cfg = [];
% cfg.fc = {'TT7_SS_01_Good.t' };
% cfg.getTTnumbers = 0;
% S = LoadSpikes(cfg);
% 
% Q = MakeQfromS([], S);
% 
% %% try a Q with Ca binaries
% 
% plot(Ca_1)
% % hold on
% % plot(ms_seg.time, ms_seg.Binary(:,1))
% 
% % Ca_Q = MakeQfromS([],Ca_TS);
% 
% % Ca_Q_SWR = MakeQfromS2([], Ca_TS, SWRs_wide)
% 
% %% get a co-occurance measure for the while block
% 
% cfg = [];
% cfg.nShuffle = 1000;
% cfg.useMask = 1;
% cfg.outputFormat = 'vectorU';
% CC = CoOccurQ(cfg,Ca_Q_SWR);
%% Prepare for SeqNMF

SWR_cat_data = [];
t_win = [0 2.5];
idx_win = ceil(t_win*(mode(diff(ms_seg.time)))); % window (in index values) around the event.

% h = imagesc(co_mat(:,:,1));

for iE =1: length(swr_ms_idx_tstart)%:-1:1 % loop SWRS
    
    if swr_ms_idx_tstart(iE)+idx_win(1) <= 0
        this_evt = ms_seg.Binary(1:swr_ms_idx_tstart(iE)+idx_win(2),:); % get all the values within this event.
        
    elseif swr_ms_idx_tstart(iE)+idx_win(2) >= length(ms_seg.Binary)
        this_evt = ms_seg.Binary(swr_ms_idx_tstart(iE)+idx_win(1):end,:); % get all the values within this event.
        
    else
        this_evt = ms_seg.Binary(swr_ms_idx_tstart(iE)+idx_win(1):swr_ms_idx_tstart(iE)+idx_win(2),:); % get all the values within this event.
    end
    SWR_activity(:,iE) = sum(this_evt,1) >0; % see if anything was active.
    
    %     corr_mat(:,:,iE) = corr(SWR_activity', 'rows', 'pairwise');
    
    for ii = length(SWR_activity(:,iE)):-1:1
        for jj = length(SWR_activity(:,iE)):-1:1
            if SWR_activity(ii,iE) == 1 && SWR_activity(jj,iE) == 1
                co_mat(ii,jj,iE) = 1;
            elseif (SWR_activity(ii,iE) == 1 && SWR_activity(jj,iE) == 0) || (SWR_activity(ii,iE) == 0 && SWR_activity(jj,iE) == 1)
                co_mat(ii,jj,iE) = 0;
            elseif SWR_activity(ii,iE) == 0 && SWR_activity(jj,iE) == 0
                co_mat(ii,jj,iE) = 0;
                
            end
        end
    end
%         h.CData = co_mat(:,:,iE);
%         drawnow
%         pause(.5)
    
    %     for iCell = length(SWR_activity(:,iE)):-1:1
    SWR_cat_data = [SWR_cat_data; this_evt];
    %     end
end


%% corrgram
  for ii = length(SWR_activity(:,iE)):-1:1

      prob_active(ii)  = nanmean(co_mat(ii,:,:), 'all');
        for jj = length(SWR_activity(:,iE)):-1:1
            corr_p_mat(ii, jj) = nanmean(co_mat(ii,jj,:), 'all');
        end
  end
  
  p0 = nanmean(SWR_cat_data,2);
  
  p1 = nan(size(SWR_activity,1));

  for ii = length(p1):-1:1
       for jj = length(p1):-1:1
       p1(ii,jj) = p0(ii).*p0(jj);
       end
  end
  
  figure(111)
  imagesc(corr_p_mat)

%% plug it into SeqNMF
addpath(PARAMS.code_seqnmf_dir)
% load('D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\Inter\537\12_5_2019_537day1\all_binary_post_SW.mat')


data_in = all_binary_post_SW'; 

data_in = SWR_cat_data'; 

data_in = ms_seg.Binary(:,1:300)';

data_in = ms_seg.RawTraces(1:1000,:)'; 

data_in = this_evt'; 
Fs =mode(diff(ms_seg.time));

%% try hat5 raw
% for ref: J:\Williams_Lab\Jisoo\Jisoo_Project\Results\PV1043\LTD5_results
load('J:\Williams_Lab\Jisoo\Jisoo_Project\RawData\pv1043\6_15_2019_PV1043_LTD5\H13_M30_S2_LTD5\ms.mat')
clearvars -except ms PARAMS
load('J:\Williams_Lab\Jisoo\Jisoo_Project\RawData\pv1043\6_15_2019_PV1043_LTD5\H13_M30_S2_LTD5\behav.mat')

ms = msExtractBinary_detrendTraces(ms);

data_in = ms.Binary(1:10000,:)'; 
pos(:,1) = interp1(behav.time,behav.position(:,1),ms.time);
pos(:,2) = interp1(behav.time,behav.position(:,2),ms.time);
velocity = interp1(behav.time, behav.speed, ms.time); 


% remove inactive cells
total_act = sum(data_in,2);

keep_idx = total_act>0;

data_in = data_in(keep_idx,:);
% limit to movement. 
% keep_idx = velocity > 1;
% 
% pos_move = pos(keep_idx);
% velo_move = velocity(keep_idx); 
% data_in = ca_binday(:,keep_idx); 
Fs = mode(diff(ms.time)); 

Ls = fliplr([0.5 1 2 10 50 100 1000]);
%%
addpath(PARAMS.code_seqnmf_dir)

cd('D:\Dropbox (Williams Lab)\Jisoo\JisooProject2020\Inter\PV1069\10_18_2019_PV1069_HATD5');

data_in = all_binary_post_REM';
data_in = all_RawTraces_post_REM';
data_in = ms_trk.Binary'; 
% % normalize(
% for ii = size(data_in,2):-1:1
%     data_out(:,ii) = (data_in(:,ii) - min(data_in(:,ii))) / ( max(data_in(:,ii))); 
% end
Fs = 33;

%% 
%% Fit with seqNMF
%% break data into training set and test set
splitN = floor(size(data_in,2)*.75); 
% splitS = floor(size(SONG,2)*.75); 
trainNEURAL = data_in(:,1:splitN); 
% trainSONG = SONG(:,1:splitS); 
testNEURAL = data_in(:,(splitN+1):end); 
% testSONG = SONG(:,(splitS+1):end); 
%% plot one example factorization
rng(235); % fixed rng seed for reproduceability
X = trainNEURAL;
for iS = length(Ls):-1:1
K = 4;
L = Ls(iS); % units of seconds
Lneural = ceil(L*Fs);  
% Lsong = ceil(L*SONGfs);
shg
display('Running seqNMF on real neural data (from songbird HVC, recorded by Emily Mackevicius, Fee Lab)')
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
all_sweeps{iS}.K = K;
all_sweeps{iS}.L = L;
all_sweeps{iS}.W = W;
all_sweeps{iS}.H = H;
all_sweeps{iS}.sig = is_significant;
all_sweeps{iS}.Train = trainNEURAL;
all_sweeps{iS}.Trest = testNEURAL;

saveas(gcf,['Seq_Sweeps' filesep 'Sweep_' num2str(L)], 'png')
close all
% clearvars -except all_sweeps trainNEURAL testNEURAL Ls

end
%%
% plot, sorting neurons by latency within each factor
[max_factor, L_sort, max_sort, hybrid] = helper.ClusterByFactor(W(:,:,:),1);
indSort = hybrid(:,3);
tstart = 180; % plot data starting at this timebin
figure; WHPlot(W(indSort,:,:),H(:,tstart:end), X(indSort,tstart:end), ...
    0,trainSONG(:,floor(tstart*Fs/VIDEOfs):end))
title('Significant seqNMF factors, with raw data')
figure; WHPlot(W(indSort,:,:),H(:,tstart:end), ...
    helper.reconstruct(W(indSort,:,:),H(:,tstart:end)),...
    0,trainSONG(:,floor(tstart*Fs/VIDEOfs):end))
title('SeqNMF reconstruction')







