% dSub cell classification sandbox

%% add paths

close all
restoredefaultpath
global PARAMS
os = computer;

if ismac
    PARAMS.data_dir = '/Users/jericcarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1'; % where to find the raw data
    PARAMS.inter_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/Users/jericcarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/Users/jericcarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    PARAMS.ft_code_dir = '/Users/jericcarmichael/Documents/GitHub/fieldtrip'; % FieldTrip toolbbox (used for spectrogram visualization)
    
elseif strcmp(os, 'GLNXA64')
    
    PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/II_classification'; % where to find the raw data
    PARAMS.inter_dir = '/home/ecarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.stats_dir = '/home/ecarmichael/Documents/Williams_Lab/Stats/'; % where to put the statistical output .txt
    PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
    PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    PARAMS.ML_spike.brick = '/home/ecarmichael/Documents/GitHub/ML_spike/brick'; % the brick part of the ML Spike spike esitmation package from Deneux et al., 2016
    PARAMS.ML_spike.spike = '/home/ecarmichael/Documents/GitHub/ML_spike/spikes'; % the spike part of the ML Spike spike esitmation package from Deneux et al., 2016

    PARAMS.seqNMF_dir = '/home/ecarmichael/Documents/GitHub/seqNMF';% seqNMF pathway. used for sequence detection.
    
else
    disp('on a PC fill this in yourself....')
end


rng(11,'twister') % for reproducibility


% add the required code
addpath(genpath(PARAMS.code_base_dir));
addpath(genpath(PARAMS.code_CEH2_dir));


% add the ML_spike
addpath(PARAMS.ML_spike.spike);
addpath(PARAMS.ML_spike.brick);

cd(PARAMS.data_dir) % move to the data folder

% try the newer NLX loaders for UNIX
[~, d] = version;
if str2double(d(end-3:end)) >2014 && strcmp(os, 'GLNXA64')
    rmpath('/Users/jericcarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared/io/neuralynx')
    addpath(genpath('/Users/jericcarmichael/Documents/NLX_loaders_UNIX_2015'))
    disp('Version is greater than 2014b on UNIX so use updated loaders found here:')
    which Nlx2MatCSC
end

clear d os


%% load the Miniscope data
this_dir = dir;

for iF = 1:length(this_dir)
   if strfind(this_dir(iF).name, 'ms.mat')
       warning('off')
       temp = load(this_dir(iF).name);
       All_sess.(['M' this_dir(iF).name(1:3)]) = temp.ms;
       clear temp
       warning('on')
       fprintf('<strong>MS_Ca_classification</strong>: %s loaded as a field in All_sess\n', this_dir(iF).name)
       
   else
       continue
   end
end


%% Get data ready for ML_spike
this_sess = All_sess.M567; % for looping later. 
%estimate framerate
dt = 1/mode(diff(this_sess.time));

%convert ms data to cells

calcium = cell(1,size(this_sess.FiltTraces,2));

for iC = size(this_sess.RawTraces,2):-1:1
    
    calcium{iC} = this_sess.RawTraces(:,iC);
    
end



%% Spike estimation (using naive parameters)
% We first apply the MLspike algorithm using some fixed parameters.
% To run the algorithm on other data, set variable 'calcium' as a vector
% time courses of the calcium signal if there is only one trial, or to a cell
% array of such time courses if there are many trials.


% parameters
% (get the default parameters)
par = tps_mlspikes('par');
% (indicate the frame duration of the data)
par.dt = dt;
% (set physiological parameters)
par.a = 0.07; % DF/F for one spike
par.tau = 1; % decay time constant (second)
par.saturation = 0.1; % OGB dye saturation
% (set noise parameters)
par.finetune.sigma = .02; % a priori level of noise (if par.finetune.sigma
                          % is left empty, MLspike has a low-level routine 
                          % to try estimating it from the data
par.drift.parameter = .01; % if par.drift parameter is not set, the 
                           % algorithm assumes that the baseline remains
                           % flat; it is also possible to tell the
                           % algorithm the value of the baseline by setting
                           % par.F0
% (do not display graph summary)
par.dographsummary = false;

ML_out.par = par;

% spike estimation
[ML_out.spikest ML_out.fit ML_out.drift] = spk_est(ML_out.calcium,ML_out.par);

figure(1)
spk_display(ML_out.dt,ML_out.spikest,{ML_out.calcium ML_out.fit ML_out.drift})
set(1,'numbertitle','off','name','MLspike alone')

%% one cell at a time?

this_Cell = 1;

figure(this_Cell)
spk_display(ML_out.dt,ML_out.spikest{this_Cell},{ML_out.calcium{this_Cell}, ML_out.fit{this_Cell}, ML_out.drift{this_Cell}})
set(1,'numbertitle','off','name','MLspike alone')

%%  use ML spike to pull out some 


load('*ML.mat')


figure(1)
spk_display(dt,spikest,{ML_out.calcium fit drift})
set(1,'numbertitle','off','name','MLspike alone')