function Beluga_gen_PPC_sl(sess_dir)

code_dir = '/home/ecar/Github/CEH2/Cluster/PPC';

if ~exist(code_dir)
    mkdir(code_dir)
end

cd(sess_dir)


if isempty(dir('*meta.m'))
    MS_Write_meta_dSub;
end

Meta = MS_Load_meta;


t_files = dir('*.t');

for iT = length(t_files):-1:1
    if str2double(t_files(iT).name(end-2)) < 5 || contains(t_files(iT).name(end-2), {'A', 'B', 'C', 'D'})    
        MS_write_PPC_slurm(code_dir, sess_dir, t_files(iT).name(1:end-2)); 

    end
end

