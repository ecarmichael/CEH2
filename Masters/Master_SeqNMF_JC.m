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
    PARAMS.raw_data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData'; % raw data location.
    PARAMS.inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\SeqNMF_EC'; % where to put intermediate files
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
% define the type of cells to use

cell_type = 'Anx';

%% set up the directories to process

session = {'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData\pv1069\10_14_2019_PV1069_HATD1',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData\pv1060\11_19_2019_PV1060_HATD1',...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData\pv1060\11_21_2019_PV1060_HATD3'}
    

% 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData\pv1069\10_18_2019_PV1069_HATD5\H13_M5_S42_HATD5',...

%     'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData\pv1069\10_22_2019_PV1069_HATSwitch\H13_M4_S44_HATD6',...
    %   'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData\pv1060\11_26_2019_PV1060_HATSwitch\H13_M5_S15_HATSwitch',...
% 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData\pv1069\7_12_2019_PV1069_LTD5\H12_M41_S29_LTD5', ...
% 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData\pv1069\7_8_2019_PV1069_LTD1\H12_M36_S0_LTD1',...
% 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData\pv1060\7_15_2019_PV1060_LTD1\H13_M22_S5_LTD1',...

% 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData\pv1060\11_19_2019_PV1060_HATD1',...
% 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData\pv1060\11_21_2019_PV1060_HATD3'};
%'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\RawData\pv1060\7_15_2019_PV1060_LTD1\H13_M22_S5_LTD1'};

pro_main = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting'; % where to find processed data

REM_main = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter';
%% cycle through some sessions

for iSess = 3:length(session)
    cfg.place = 1;
    cfg.anx = 0;
    MS_Run_REM_SeqNMF(cfg, session{iSess}, pro_main, REM_main)
    close all
    
%     if contains(session{iSess}, 'HAT')
%         cfg.place = 0;
%         cfg.anx = 1;
%         MS_Run_REM_SeqNMF(cfg, session{iSess}, pro_main, REM_main)
%         close all
%     end
end