
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
    if strcmpi(getenv('USERNAME'), 'ecarmichael')
        PARAMS.data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/JC/7_12_2019_PV1069_LTD5'; % where to find the raw data
        PARAMS.code_base_dir = '/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'; % where the codebase repo can be found
        PARAMS.code_CEH2_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % where the multisite repo can be found
    elseif strcmpi(getenv('USERNAME'), 'williamslab')
        PARAMS.data_dir = '/home/williamslab/Dropbox (Williams Lab)/10.Manifold/pv1043/LTD1'; % where to find the raw data
        PARAMS.code_CEH2_dir = '/home/williamslab/Documents/Github/CEH2'; % where the multisite repo can be found
        PARAMS.code_seqNMF = '/home/williamslab/Documents/Github/seqNMF'; 
    end
else
    PARAMS.data_dir = 'D:\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\Processed place'; % where to find the raw data
    PARAMS.inter_dir = 'D:\Dropbox (Williams Lab)\Jisoo\Upload\pv1002\Inter'; % where to put intermediate files
    PARAMS.code_CEH2_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2'; % where the multisite repo can be found
    PARAMS.code_seqnmf_dir = 'C:\Users\ecarm\Documents\GitHub\seqNMF'; % where the multisite repo can be found
end

rng(11,'twister') % for reproducibility

% add the required code
addpath(genpath(PARAMS.code_CEH2_dir));
cd(PARAMS.data_dir) % move to the data folder

clear d os

%% load some sample data
load('all_binary_post_REM.mat')
load('ms.mat')
load('behav.mat')
load('decoding.mat')

%% prepare for Seq

Fs =mode(diff(ms.time));

this_data = all_binary_post_REM;
this_pos = decoding.REM_decoded_probabilities;

PosFs = Fs; % b/c we already interpolated the position data to the ms.time. 

% break data into training set and test set
splitN = floor(size(this_data,2)*.75); 
splitS = floor(size(this_pos,2)*.75);

% make the test/trian sets.  
trainNEURAL = this_data(:,1:splitN); 
trainPOS = this_pos(:,1:splitS); 
testNEURAL = this_data(:,(splitN+1):end); 
testPOS = this_pos(:,(splitS+1):end); 
%% plot one example factorization
rng(235); % fixed rng seed for reproduceability
X = trainNEURAL;
figure
% loop L values
% Ls = fliplr([0.5 1 2 10 50 100 1000]);
% 
% for iS = length(Ls):-1:1
K = 2;
L = 10;
% L = Ls(iS); % units of seconds
Lneural = ceil(L*Fs);  
% Lsong = ceil(L*SONGfs);
shg
subplot(2,2,1)
display('Running seqNMF...')
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





