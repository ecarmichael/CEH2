function Beluga_PPC(sess_dir, t_name)

addpath(genpath('/home/ecar/Github/CEH2'));
addpath(genpath('/home/ecar/Github/vandermeerlab')); 
addpath('/home/ecar/Github/fieldtrip')
ft_defaults


disp('Codebases loaded')

data_dir = '/lustre04/scratch/ecar';

cd([data_dir filesep sess_dir(1:3) filesep sess_dir])

if isempty(dir('*meta.m'))
    MS_Write_meta_dSub;
end

Meta = MS_Load_meta;
%% get the list of good cells

   MS_get_PPC([],t_name, {Meta.goodCSC}); 
       
       disp('Complete')
        
