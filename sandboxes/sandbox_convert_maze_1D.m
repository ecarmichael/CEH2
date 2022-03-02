%% sandbox_convert MAze to 1D

addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'));
addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));

cd('C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3');

load('behav_DLC.mat');
load('ms_trk.mat'); 

%% office computer

addpath(genpath('/home/williamslab/Documents/Github/CEH2'));
addpath(genpath('/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared'));

cd('/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/Maze_Ca/inter/pv1254/2021_12_17_pv1254_MZD3');

load('behav_DLC.mat');
load('ms_trk.mat'); 
%% play with linerization using one direction at a time
pos = MS_behav2tsd(behav); 
pos.tvec = pos.tvec/1000; 
% limit to movement
linspeed = getLinSpd([],pos); % linear speed

linspeed.data = smooth(linspeed.data, floor(1/mode(diff(linspeed.tvec)))*.5);
move_idx = linspeed.data > 3; 

SCoorD_L = StandardizeCoord([], behav.CoorD_L, 100);
SCoorD_R = StandardizeCoord([], behav.CoorD_R, 100);


