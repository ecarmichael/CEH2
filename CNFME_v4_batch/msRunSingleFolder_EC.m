data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Morph\ECQL\ECQL\Experiment0\dSub_g8_B\customEntValHere\2023_03_08\14_50_30\Grey';

% numcores = feature('numcores') 

% Create a "local" cluster object
local_cluster = parcluster('local');

% Modify the JobStorageLocation to $SLURM_TMPDIR
% local_cluster.JobStorageLocation = getenv(cd)

% Start the parallel pool
parpool(local_cluster,2);

analysis_path = data_dir; %put in path to folder in question here (or make copy of this file within the folder you want to run and enjoy)

addpath(genpath(['C:\Users\ecarm\Documents\GitHub\Brandon-Williams-Analysis-Package\']))
rmpath(genpath(replace(['C:\Users\ecarm\Documents\GitHub\Brandon-Williams-Analysis-Package/Imported Analysis Scripts/AlexAnalysis/'], {'\', '/'}, {filesep, filesep})));
rmpath(genpath(replace(['C:\Users\ecarm\Documents\GitHub\Brandon-Williams-Analysis-Package/Imported Analysis Scripts/CNMF_E-master/'], {'\', '/'}, {filesep, filesep})));

oldcd = analysis_path;

cd(replace('C:\Users\ecarm\Documents\GitHub\Brandon-Williams-Analysis-Package/Imported Analysis Scripts/CNMF_E-master_Linux/', {'\', '/'}, {filesep, filesep}))

cvx_setup

cd(oldcd)

msRun2020_newSoft_EC(analysis_path)
