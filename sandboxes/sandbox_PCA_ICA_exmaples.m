%% sandbox_PCA_ICA_examples



%% %%%%% EPHYS %%%%%%%%%%%%
load('/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Radial/Inter/M30_2023_06_10_Rad5.mat')


%% get some assembly patterns

[A_temp, T_proj] = MS_PCA_ICA([], out.Encode.S, out.Encode.pos, out.Encode.csc)