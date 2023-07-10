%%sandbox_Master_NOL

kilo_dir = '/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Eric/NOL/M30_2023-06-27_NOL1/kilo_out';

OE_dir = '/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Eric/NOL/M30_2023-06-27_NOL1/Record Node 120'; 

DLC_dir = {'/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Eric/NOL/M30_2023-06-27_NOL1/2023_06_27/09_24_45/My_WebCam'; ...
    '/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Eric/NOL/M30_2023-06-27_NOL1/2023_06_27/14_14_54/My_WebCam'}; 

save_dir = '/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Eric/NOL/inter'; 


if strcmp(computer, 'GLNXA64')
    
    addpath(genpath('/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared'));
    
    addpath(genpath('/home/williamslab/Documents/Github/CEH2'))
    
    addpath(genpath('/home/williamslab/Documents/Github/npy-matlab'))

end





%%
Master_NOL_preprocess([], kilo_dir, OE_dir, DLC_dir, save_dir)





%%  REM rad

data_dir = '/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Radial/Radial_ephys/M30/M30_2023-06-06_Rad1'; 


if strcmp(computer, 'GLNXA64')
    
    addpath(genpath('/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared'));
    
    addpath(genpath('/home/williamslab/Documents/Github/CEH2'))
    
    addpath(genpath('/home/williamslab/Documents/Github/npy-matlab'))

end

%% 

out = Master_REM_radial_preprocess([], data_dir)