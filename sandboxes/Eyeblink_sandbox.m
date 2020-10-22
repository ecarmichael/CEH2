%% Eyeblink Gamma Spec sandbox


close all
restoredefaultpath
global PARAMS

os = computer;

if ismac
    PARAMS.data_dir = '/Volumes/Fenrir/State_dep'; % where to find the raw data
    PARAMS.inter_dir = '/Volumes/Fenrir/State_dep/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/Volumes/Fenrir/State_dep/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_state_dir = '/Users/jericcarmichael/Documents/GitHub/EC_State'; % where the multisite repo can be found
    PARAMS.ft_dir = '/Users/jericcarmichael/Documents/GitHub/fieldtrip'; % if needed. 
    PARAMS.filter_dir = '/Volumes/Fenrir/State_dep/Filters/'; % stores custom built filters for speed.  
    
elseif strcmp(os, 'GLNXA64')
    
    PARAMS.data_dir = '/media/ecarmichael/Fenrir/Eyeblink/M142-2020-10-05-CDOD6'; % where to find the project data
    PARAMS.raw_data_dir = '/media/ecarmichael/Fenrir/Eyeblink/M142-2020-10-05-CDOD6'; % where to find the raw data
    PARAMS.inter_dir = '/media/ecarmichael/Fenrir/Eyeblink/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/media/ecarmichael/Fenrir/Eyeblink/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_state_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    PARAMS.ft_dir = '/home/ecarmichael/Documents/GitHub/fieldtrip'; % if needed. 
    PARAMS.filter_dir = '/media/ecarmichael/Fenrir/State_dep/Filters/'; % stores custom built filters for speed.
    
else
    disp('on a PC')
%     PARAMS.data_dir = 'G:\Multisite\'; % where to find the raw data
%     PARAMS.inter_dir = 'G:\Multisite\temp\'; % where to put intermediate files
%     PARAMS.stats_dir = 'G:\Multisite\temp\Stats\'; % where to put the statistical output .txt
%     PARAMS.code_base_dir = 'D:\Users\mvdmlab\My_Documents\GitHub\vandermeerlab\code-matlab\shared'; % where the codebase repo can be found
%     PARAMS.code_MS_dir = 'D:\Users\mvdmlab\My_Documents\GitHub\EC_Multisite'; % where the multisite repo can be found
end
rng(11,'twister') % for reproducibility


% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_state_dir));
cd(PARAMS.data_dir) % move to the data folder


%% Load data

cfg = [];
cfg.fc = {'CSC40.ncs'};
cfg.desired_sampling_frequency = 2000; % decimate to 2k
csc = MS_LoadCSC(cfg);

%% load events

evts = LoadEvents([])

%% plot the events in time to check when things happen
c_ord = linspecer(length(evts.t)); 

figure(1)
plot(csc.tvec - csc.tvec(1), NaN(1,length(csc.tvec)))

hold on
line([evts.t{1}; evts.t{1}], [ones(1,length(evts.t{1})); 1.3*ones(1,length(evts.t{1}))], 'color', c_ord(2,:))
line([evts.t{2}; evts.t{2}], [ones(1,length(evts.t{2})); 1.3*ones(1,length(evts.t{2}))], 'color', 'k')

% TTL events
labels{1} = 'Start';
labels{2} = 'Stop'; 
for iE = 3:length(evts.t)
    plot(evts.t{iE}, 1+(iE*0.02)*ones(1,length(evts.t{iE})), '*','color', c_ord(iE,:));
    labels{iE} = num2str(iE); 
end

legend(labels)
%% Restrict to whatever the third recording block is.  It seems to be the same length as the 4th.  Maybe CS+/CS- pairs swap?


csc_r = restrict(csc, evts.t{1}(3), evts.t{2}(6));

evts_r = restrict(evts, evts.t{1}(3), evts.t{2}(6)); 

%% extract the CS+ and CS- trials based on the context from Jimmie
% 08 is rewarded when preceded by 04.
% 02 is rewarded when preceded by 40
% A trial goes, 1 s context cue (either 04 or 40), 2 s delay, 1 s outcome-predictive cue (08 or 02)"

% cont_A = 7; rew_A = 8;
% cont_B = 13; rew_B = 6; 

cont_A = evts_r.t{7}; 
cont_B = evts_r.t{13};

CS_A = evts_r.t{8};
CS_B = evts_r.t{6};

Rew = evts_r.t{12}; % rewards
Air = evts_r.t{9}; % air puffs
Licks = evts_r.t{4}; %licks
                                                              
Rew_AA = []; Rew_AB = [];
Miss_AA =[]; Miss_AB = [];
Lick_AA =[]; Lick_AB = [];
Con_AA = []; Con_AB = [];
%%RETHINK THIS YOU WANT THE CS TIMES NOT THE CON TIMES!!!! ADD IN REWARDED
%%VS MISS< THESE ONLY CONTEXT+STIMULUS NOT ACTUAL REWARD

% extract trials for context A
for iA = 1:length(CS_A)
    all_CS_AA(iA) = cont_A(nearest_idx3(CS_A(iA),cont_A)) - CS_A(iA);
    prox_lick(iA) = Licks(nearest_idx3(CS_A(iA),Licks)) - CS_A(iA);
    if abs(all_CS_AA(iA)) <= 4
        Con_AA = [Con_AA, all_CS_AA(iA)];
        Rew_AA = [Rew_AA, CS_A(iA)];
        if prox_lick(iA) >0  && prox_lick(iA) < 2
            Lick_AA = [Lick_AA, CS_A(iA)];
        end
    end
end
% extract context A CS B
for iA = 1:length(CS_A)
    all_CS_AB(iA) = cont_B(nearest_idx3(CS_A(iA),cont_B)) - CS_A(iA);
    prox_lick(iA) = Licks(nearest_idx3(CS_A(iA),Licks)) - CS_A(iA);
    if abs(all_CS_AB(iA)) <= 4
        Con_AB = [Con_AB, all_CS_AB(iA)];
        Rew_AB = [Rew_AB, CS_A(iA)];
        if prox_lick(iA) >0  && prox_lick(iA) < 2
            Lick_AB = [Lick_AB, CS_A(iA)];
        end
    end
end

Rew_BA = []; Rew_BB = [];
Miss_BA =[]; Miss_BB = [];
Lick_BA =[]; Lick_BB = [];
Con_BA = []; Con_BB = [];

% extract trials for context B CS A
for iA = 1:length(CS_A)
    all_CS_BA(iA) = cont_B(nearest_idx3(CS_A(iA),cont_B)) - CS_A(iA);
    prox_lick(iA) = Licks(nearest_idx3(CS_A(iA),Licks)) - CS_A(iA);
    if abs(all_CS_BA(iA)) <= 4
        Con_BA = [Con_BA, all_CS_BA(iA)];
        Rew_BA = [Rew_BA, CS_A(iA)];
        if prox_lick(iA) >0  && prox_lick(iA) < 2
            Lick_BA = [Lick_BA, CS_A(iA)];
        end
    end
end
% extract context B CS B
for iA = 1:length(CS_B)
    all_CS_BB(iA) = cont_B(nearest_idx3(CS_B(iA),cont_B)) - CS_B(iA);
    prox_lick(iA) = Licks(nearest_idx3(CS_B(iA),Licks)) - CS_B(iA);
    if abs(all_CS_BB(iA)) <= 4
        Con_BB = [Con_BB, all_CS_BB(iA)];
        Rew_BB = [Rew_BB, CS_B(iA)];
        if prox_lick(iA) >0  && prox_lick(iA) < 2
            Lick_BB = [Lick_BB, CS_B(iA)];
        end
    end
end
%% prepare the data for FT

addpath(PARAMS.ft_dir)
ft_defaults


subplot(3,2,1)
Triggered_Spec_FT(csc_r, Rew_AA, 'AA Reward')

subplot(3,2,2)
Triggered_Spec_FT(csc_r, Rew_BB, 'BB Reward')

subplot(3,2,3)
% Triggered_Spec_FT(csc_r, Lick_AA, 'Lick AA')
Triggered_Spec_FT(csc_r, Con_AA, 'Lick AA')

subplot(3,2,4)
% Triggered_Spec_FT(csc_r, Lick_AB, 'Lick AB')
Triggered_Spec_FT(csc_r, Con_AB, 'Lick AB')

subplot(3,2,5)
% Triggered_Spec_FT(csc_r, Lick_BA, 'Lick BA')
Triggered_Spec_FT(csc_r, Con_BA, 'Lick BA')

subplot(3,2,6)
% Triggered_Spec_FT(csc_r, Lick_BB, 'Lick BB')
Triggered_Spec_FT(csc_r, Con_BB, 'Lick BB')

%%
loop_n  = 0;
% for iE = [8, 6, 9, 12]
%     loop_n = loop_n+1; 
subplot(2,2,1)
Triggered_Spec_FT(csc_r, Rew_AA, 'AA_correct')
% end


