%% CA2 + ephys sandbox

%% add paths

close all
restoredefaultpath
global PARAMS

if isunix
    PARAMS.data_dir = '/Users/jericcarmichael/Documents/Williams_Lab/7_12_2019_PV1069_LTD5'; % where to find the raw data
    PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/Users/jericcarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
else
 disp('on a PC')
end


rng(11,'twister') % for reproducibility


% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_CEH2_dir));
cd(PARAMS.data_dir) % move to the data folder


%%  Laod some stuff

load('ms.mat');
load('SFP.mat');


%% make a video

for iframe =  1:size(SFP, 3)
    
    imagesc(SFP(:,:,iframe))
    M(iframe) = getframe;
end


%% GET nlx data
cfg = [];
cfg.fc = {'CSC1.ncs', 'CSC8.ncs'};
csc = LoadCSC(cfg);

cfg = [];
evt = LoadEvents(cfg);

%restrict data
csc_r = restrict(csc, evt.t{3}(1), evt.t{3}(end));
%% plot
figure(8)

% plot(csc.tvec(1:10000), csc.data(1,1:10000), csc.tvec(1:10000), csc.data(2,1:10000),evt.t{3}, '*k' )
t_start = nearest_idx3(csc.tvec, evt.t{3}(1));
t_end = nearest_idx3(csc.tvec, evt.t{3}(end));

%%
plot(csc_r.tvec, csc_r.data(1,:),...
    csc_r.tvec,csc_r.data(2,:))%,...

hold on
plot(evt.t{3},max(csc_r.data(1,:)), '*k' )
plot(evt.t{4},max(csc_r.data(1,:)), '*c' )


