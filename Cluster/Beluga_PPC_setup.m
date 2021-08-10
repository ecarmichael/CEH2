function Beluga_PPC_setup(sess_dir)

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

s_files = dir('*.t');

for iS = length(s_files):-1:1
    if str2double(s_files(iS).name(end-2)) < 5 || contains(s_files(iS).name(end-2), {'A', 'B', 'C', 'D'})

       MS_get_PPC([], s_files(iS).name, Meta.goodCSC); 
        
    end
end



