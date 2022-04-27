%% seq NMF JC test
clear all; close all; 
addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared')); % where the codebase repo can be found
addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2')); % where the multisite repo can be found
addpath(genpath('C:\Users\ecarm\Documents\GitHub\seqNMF'));
data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1060\HATD1';

cd(data_dir)

load('behav.mat');
load('ms.mat', 'ms')

%% try it with the CC seqNMF script.

k = 4;
l = 2;

close all
CC_SeqNMF(ms, behav, k, l, 0.0034)



%% same but REM
l = 2;

load('all_binary_post_REM.mat')

CC_SeqNMF_Sleep(all_binary_post_REM,33, k, l, 0.0034)



%% comapre the Seq templates with the detected replays
load('C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold\pv1060\HATD1\Seq_out_sleep_L0p5\Seq_out_sleep.mat'); 

load('C:\Users\ecarm\Dropbox (Williams Lab)\Decoding\pv1060\HATD1\1000shuffling.mat','Replay_score_actual', 'Final_Replay_score', 'Final_start_frame');
% replay_evts = load(

%% try a peri event plot of the H mat values from Seq
 iter = 1;
 this_W = Seq_out_sleep.W{iter}; 
 this_H = Seq_out_sleep.H{iter};

 % test figure
 c_ord = linspecer(size(this_H,1)); 
 figure(101)
 clf
 subplot(5,1,1:2)
  hold on

 for ii = size(this_H,1):-1:1
     plot(1:length(this_H(ii,:)), (this_H(ii,:)./max(this_H(ii,:)))-ii, 'color', c_ord(ii,:));  
 end
xlim([1 size(this_H,2)])


%% loop for SCEs
data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold';
cd(data_dir)
subs = dir('pv*');


for iSub = 1:length(subs)
    if strcmp(subs(iSub).name, 'pv1215')
        continue
    end
    cd([data_dir filesep subs(iSub).name])
    
    sess = {'LTD1','LTD3','LTD5','HATD1','HATD3','HATD5','HATDSwitch'}; 

    
    for iSess = 1:length(sess)
        if ~exist([data_dir filesep subs(iSub).name filesep sess{iSess}], 'dir')
            continue
        else
            cd([data_dir filesep subs(iSub).name filesep sess{iSess}])
            if ~exist('behav.mat')
                continue
            end
            load('behav.mat')
            load('ms.mat')
            
            %             All.(subs(iSub).name).(sess(iSess).name) = MS_get_SCE_JC(ms, behav);
            All_out{iSub,iSess} = MS_get_SCE_JC(ms, behav);
            
            
            clear behav ms
        end
    end
    
    %         if isfield(All.(subs(iSub).name), 'LTD1')
    %         LTDs = fieldnames
    %
    %
    %
    %         all_LTD(iSub) = All.(subs(iSub).name).LTD1
    
end

%% collect and plot

for iSess = 1:length(sess)
    this_bar = [];
    flag = []; 
    for iSub = 1:length(subs)
        if isempty(All_out{iSub, iSess})
%             this_bar{iSub} = [];
            flag(iSub) = 1; 
            continue
        else
            flag(iSub) = 0; 
            this_bar(iSub,:) = All_out{iSub, iSess}.rate_vec;
            bins = All_out{iSub, iSess}.bins;
        end
    end
    
    this_bar(logical(flag),:) = []; 
    All_rates.(sess{iSess}) = this_bar;
end

%%
cord = linspecer(20);
figure(101)
clf
m = 2; n = 4;
subplot(m, n, 1)
hold on
h =bar(bins, nanmean(All_rates.LTD1,1), 'facecolor',cord(1,:));
errorbar(bins, nanmean(All_rates.LTD1,1), nanstd(All_rates.LTD1,1),'k', 'LineStyle', 'none');
set(h,'BarWidth',1,'EdgeColor','none');
title(['LTD1 (n:' num2str(size(All_rates.LTD1,1)) ')']); 
ylabel('SCE rate /s')
ylim([0 max(nanmean(All_rates.LTD1,1)+nanstd(All_rates.LTD1,1))])
y_tick = get(gca, 'ytick');
rectangle('Position', [0, max(nanmean(All_rates.LTD1,1)+nanstd(All_rates.LTD1,1))*.9,100,max(nanmean(All_rates.LTD1,1)+nanstd(All_rates.LTD1,1))*.1], 'FaceColor', [.5 .5 .5 .2], 'edgeColor', [.5 .5 .5 .2]);
xlim([0 100])

subplot(m, n, 2)
hold on
h=bar(bins, nanmean(All_rates.LTD3,1), 'facecolor',cord(2,:));
errorbar(bins, nanmean(All_rates.LTD3,1), nanstd(All_rates.LTD3,1),'k', 'LineStyle', 'none')
set(h,'BarWidth',1,'EdgeColor','none')
title('LTD3'); 
ylim([0 max(nanmean(All_rates.LTD3,1)+nanstd(All_rates.LTD3,1))])
y_tick = get(gca, 'ytick');
rectangle('Position', [0, max(nanmean(All_rates.LTD3,1)+nanstd(All_rates.LTD3,1))*.9,100,max(nanmean(All_rates.LTD3,1)+nanstd(All_rates.LTD3,1))*.1], 'FaceColor', [.5 .5 .5 .2], 'edgeColor', [.5 .5 .5 .2]);
xlim([0 100])

subplot(m, n, 3)
hold on
h =  bar(bins, nanmean(All_rates.LTD5,1), 'facecolor',cord(3,:));
errorbar(bins, nanmean(All_rates.LTD5,1), nanstd(All_rates.LTD5,1),'k', 'LineStyle', 'none')
set(h,'BarWidth',1,'EdgeColor','none')
title('LTD5'); 
ylim([0 max(nanmean(All_rates.LTD5,1)+nanstd(All_rates.LTD5,1))])
y_tick = get(gca, 'ytick');
rectangle('Position', [0, max(nanmean(All_rates.LTD5,1)+nanstd(All_rates.LTD5,1))*.9,100,max(nanmean(All_rates.LTD5,1)+nanstd(All_rates.LTD5,1))*.1], 'FaceColor', [.5 .5 .5 .2], 'edgeColor', [.5 .5 .5 .2]);
xlim([0 100])

subplot(m, n, 5)
hold on
h =  bar(bins, nanmean(All_rates.HATD1,1),'facecolor',cord(end-2,:));
errorbar(bins, nanmean(All_rates.HATD1,1), nanstd(All_rates.HATD1,1),'k', 'LineStyle', 'none')
set(h,'BarWidth',1,'EdgeColor','none')
title('HATD1'); 
ylabel('SCE rate /s')
xlabel('position on track (cm)')
ylim([0 max(nanmean(All_rates.HATD1,1)+nanstd(All_rates.HATD1,1))])
y_tick = get(gca, 'ytick');
rectangle('Position', [50, max(nanmean(All_rates.HATD1,1)+nanstd(All_rates.HATD1,1))*.9,50,max(nanmean(All_rates.HATD1,1)+nanstd(All_rates.HATD1,1))*.1], 'FaceColor', [.5 .5 .5 .2], 'edgeColor', [.5 .5 .5 .2]);
xlim([0 100])

subplot(m, n, 6)
hold on
h =  bar(bins, nanmean(All_rates.HATD3,1),'facecolor',cord(end-1,:));
errorbar(bins, nanmean(All_rates.HATD3,1), nanstd(All_rates.HATD3,1),'k', 'LineStyle', 'none')
set(h,'BarWidth',1,'EdgeColor','none')
title('HATD3'); 
ylim([0 max(nanmean(All_rates.HATD3,1)+nanstd(All_rates.HATD3,1))])
y_tick = get(gca, 'ytick');
rectangle('Position', [50, max(nanmean(All_rates.HATD3,1)+nanstd(All_rates.HATD3,1))*.9,50,max(nanmean(All_rates.HATD3,1)+nanstd(All_rates.HATD3,1))*.1], 'FaceColor', [.5 .5 .5 .2], 'edgeColor', [.5 .5 .5 .2]);
xlim([0 100])

subplot(m, n, 7)
hold on
h = bar(bins, nanmean(All_rates.HATD5,1),'facecolor',cord(end,:));
errorbar(bins, nanmean(All_rates.HATD5,1), nanstd(All_rates.HATD5,1),'k', 'LineStyle', 'none')
set(h,'BarWidth',1,'EdgeColor','none')
title('HATD5'); 
ylim([0 max(nanmean(All_rates.HATD5,1)+nanstd(All_rates.HATD5,1))])
y_tick = get(gca, 'ytick');
rectangle('Position', [50, max(nanmean(All_rates.HATD5,1)+nanstd(All_rates.HATD5,1))*.9,50,max(nanmean(All_rates.HATD5,1)+nanstd(All_rates.HATD5,1))*.1], 'FaceColor', [.5 .5 .5 .2], 'edgeColor', [.5 .5 .5 .2]);
xlim([0 100])

subplot(m, n, 8)
hold on
h = bar(bins, nanmean(All_rates.HATDSwitch,1),'facecolor',cord(end-6,:));
errorbar(bins, nanmean(All_rates.HATDSwitch,1), nanstd(All_rates.HATDSwitch,1),'k', 'LineStyle', 'none')
set(h,'BarWidth',1,'EdgeColor','none')
title('HATDS'); 
ylim([0 max(nanmean(All_rates.HATDSwitch,1)+nanstd(All_rates.HATDSwitch,1))])
y_tick = get(gca, 'ytick');
rectangle('Position', [0, max(nanmean(All_rates.HATDSwitch,1)+nanstd(All_rates.HATDSwitch,1))*.9,50,max(nanmean(All_rates.HATDSwitch,1)+nanstd(All_rates.HATDSwitch,1))*.1], 'FaceColor', [.5 .5 .5 .2], 'edgeColor', [.5 .5 .5 .2]);
xlim([0 100])


