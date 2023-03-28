#!/bin/bash
#SBATCH --job-name=TASKNAME
#SBATCH --account=rpp-markpb68
#SBATCH --time=8:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=14
#SBATCH --mem=285000
#SBATCH --mail-user=MYEMAIL
#SBATCH --mail-type=ALL

for i in *.avi;                                                                                                                                                                                                                                
  do name=`echo "$i" | cut -d'.' -f1`                                                                                                                                                                                                          
  echo "$name"                                                                                                                                                                                                                                
  ffmpeg -i "$i" -vcodec rawvideo "${name}_grey.avi"                                                                                                                                                                                          
  rm "${name}.avi"                                                                                                                                                                                                                          
done                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

suffix="_grey"                                                                                                                                                                                                                              
for i in *.avi;                                                                                                                                                                                                                                
  do name=`echo "$i" | cut -d'.' -f1`                                                                                                                                                                                                          
  echo "$name"                                                                                                                                                                                                                                
  cname="${name}"                                                                                                                                                                                                                              
  nname=${cname%$suffix}                                                                                                                                                                                                                      
  cp "${name}.avi" "${nname}.avi"                                                                                                                                                                                                              
  rm "${name}.avi"                                                                                                                                                                                                                          
done

sleep 10

module load StdEnv/2020
module load gentoo/2020
module load matlab/2020a

matlab -nodisplay -batch "msRunSingleFolder_Quinn"