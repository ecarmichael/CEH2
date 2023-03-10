numcores = feature('numcores') 

% Create a "local" cluster object
local_cluster = parcluster('local')

% Modify the JobStorageLocation to $SLURM_TMPDIR
local_cluster.JobStorageLocation = getenv('SLURM_TMPDIR')

% Start the parallel pool
parpool(local_cluster,numcores);

analysis_path = pwd; %put in path to folder in question here (or make copy of this file within the folder you want to run and enjoy)

addpath(genpath(['/lustre03/project/6049321/m3group/Wilson/Brandon-Williams-Analysis-Package/']))
rmpath(genpath(['/lustre03/project/6049321/m3group/Wilson/Brandon-Williams-Analysis-Package/Imported Analysis Scripts/AlexAnalysis/']));
rmpath(genpath(['/lustre03/project/6049321/m3group/Wilson/Brandon-Williams-Analysis-Package/Imported Analysis Scripts/CNMF_E-master/']));

oldcd = analysis_path;

cd('/lustre03/project/6049321/m3group/Wilson/Brandon-Williams-Analysis-Package/Imported Analysis Scripts/CNMF_E-master_Linux/')

cvx_setup

cd(oldcd)

msRun2020_newSoft_Quinn(analysis_path)
