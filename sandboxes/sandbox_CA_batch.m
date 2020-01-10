%% Ca2 processing sandbox
% ATM just used for identifying paramters and desired changes.  

%% set up 


%% add paths

close all
restoredefaultpath
global PARAMS
os = computer; 

if ismac
 error('on a PC...')

elseif strcmp(os, 'GLNXA64')

    PARAMS.raw_data_dir = '/home/ecarmichael/Documents/Williams_Lab/Raw_data/'; % where to find the raw data
    PARAMS.inter_dir = '/home/ecarmichael/Documents/Williams_Lab/Temp/'; % where to put intermediate files
    PARAMS.GE_code_dir = '/home/ecarmichael/Documents/GitHub/MiniscopeAnalysis'; % Guillaume Etter's processing codebase with core pipeline. https://github.com/etterguillaume/MiniscopeAnalysis.git
    PARAMS.cnfm_code_dir = '/home/ecarmichael/Documents/GitHub/CNMF_E'; % CNFM_E code https://github.com/zhoupc/CNMF_E
    PARAMS.CellReg_code_dir = '/home/ecarmichael/Documents/GitHub/CellReg'; % CellReg codebase. https://github.com/zivlab/CellReg
    PARAMS.NORMCOR_code_dir = '/home/ecarmichael/Documents/GitHub/NoRMCorre'; % NoRMCorre codebase. https://github.com/flatironinstitute/NoRMCorre
    PARAMS.CEH2_code_dir = '/home/ecarmichael/Documents/GitHub/CEH2'; % EC_calcium codebase. 


else
 error('on a PC...')
end


rng(11,'twister') % for reproducibility

% add the required code
addpath(genpath(PARAMS.CEH2_code_dir));
addpath(genpath(PARAMS.cnfm_code_dir));
addpath(genpath(PARAMS.GE_code_dir));
addpath(genpath(PARAMS.CellReg_code_dir));
addpath(genpath(PARAMS.NORMCOR_code_dir));

cnmfe_setup

cd(PARAMS.raw_data_dir) % move to the data folder


%% inital stuff

% temporary just use a single data set
to_process = dir(PARAMS.raw_data_dir);


all_data_dir = {};
for iF = 1:length(to_process) 
    if strcmp(to_process(iF).name, '.') || strcmp(to_process(iF).name, '..')
        continue
    else
        all_data_dir = cat(1, all_data_dir, to_process(iF).name);
    end
end


for iSess = all_data_dir
    fprintf('\nMS_Batch_Miniscope:  Processing session: %s.....\n', iSess{1})
    
    % move to the session dir
    cd([PARAMS.raw_data_dir  filesep  iSess{1}])
    
    
    
    
    
    
    
    
    
    
    
    
    
end





