function MS_write_PPC_slurm(sl_dir, sess_dir, t_name)

parts = strsplit(sess_dir, filesep);
sess_name = parts{end}; 

cd(sl_dir); 

fid = fopen(['PPC_' strrep(sess_name, '-', '_') '_' strrep(t_name, '-', '_') '.sl'], 'w+'); 


fprintf(fid, '#!/bin/bash\n#SBATCH --job-name=PPC_%s_%s', strrep(sess_name, '-', '_'), strrep(t_name, '-', '_')); 
fprintf(fid, '\n#SBATCH --account=def-wilsyl\n');
fprintf(fid, '#SBATCH --time=12:00:00\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=8\n#SBATCH --mem=32000');
fprintf(fid, '#SBATCH --mail-user=ecarmichael@gmail.com\n#SBATCH --mail-type=ALL');

fprintf(fid,'\nmodule load StdEnv/2020\nmodule load matlab/2020a\n');

fprintf(fid, '\nmatlab -nodisplay -batch "Beluga_PPC(''%s'', ''%s'')\n', sess_dir, t_name); 
fclose(fid);
