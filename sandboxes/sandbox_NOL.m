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

save_dir = '/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Radial/Inter'; 


if strcmp(computer, 'GLNXA64')
    
    addpath(genpath('/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared'));
    
    addpath(genpath('/home/williamslab/Documents/Github/CEH2'))
    
    addpath(genpath('/home/williamslab/Documents/Github/npy-matlab'))

end

%% 

out = Master_REM_radial_preprocess([], data_dir)

save([save_dir filesep strrep([out.meta.subject '_' out.meta.date '_' out.meta.session], '-', '_') '.mat'],'out');  

%% try it for a session


maps = MS_get_place_field([], restrict(out.Encode.S, out.Encode.trials(1,:), out.Encode.trials(2,:)), restrict(out.Encode.pos, out.Encode.trials(1,:), out.Encode.trials(2,:)), out.Encode.speed)


%%
maps = MS_get_place_field([], restrict(out.Recall.S, out.Recall.trials(1,:), out.Recall.trials(2,:)), restrict(out.Recall.pos, out.Recall.trials(1,:), out.Recall.trials(2,:)), out.Recall.speed)

