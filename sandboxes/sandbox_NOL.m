%%sandbox_Master_NOL

if strcmp(computer, 'GLNXA64')
    
    addpath(genpath('/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared'));
    
    addpath(genpath('/home/williamslab/Documents/Github/CEH2'))
    
    addpath(genpath('/home/williamslab/Documents/Github/npy-matlab'))
    
    pre_fix = '/home/williamslab/';
else
    pre_fix = 'C:\Users\ecarm\'; 
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'));
    
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'))
    
    addpath(genpath('C:\Users\ecarm\Documents\GitHub\npy-matlab'))
end

kilo_dir = [pre_fix strrep('Williams Lab Dropbox/Williams Lab Team Folder/Eric/NOL/M30_2023-06-27_NOL1/kilo_out', '/', filesep)];

OE_dir = [pre_fix strrep('Williams Lab Dropbox/Williams Lab Team Folder/Eric/NOL/M30_2023-06-27_NOL1/Record Node 120', '/', filesep)]; 

DLC_dir = {[pre_fix strrep('Williams Lab Dropbox/Williams Lab Team Folder/Eric/NOL/M30_2023-06-27_NOL1/2023_06_27/09_24_45/My_WebCam', '/', filesep)]; ...
    [pre_fix strrep('Williams Lab Dropbox/Williams Lab Team Folder/Eric/NOL/M30_2023-06-27_NOL1/2023_06_27/14_14_54/My_WebCam', '/', filesep)]}; 

save_dir = [pre_fix strrep('Williams Lab Dropbox/Williams Lab Team Folder/Eric/NOL/inter', '/', filesep)];




%%
out = Master_NOL_preprocess([], kilo_dir, OE_dir, DLC_dir, save_dir)




save([save_dir filesep strrep([out.meta.subject '_' out.meta.date '_' out.meta.session], '-', '_') '.mat'],'out');  


%% Raster example
enc.S = restrict(out.S, out.csc.tvec(1), out.csc.tvec(2000*10));

cfg_in = []; 
cfg_in.spkColor  = parula(length(enc.S.t)+20);
MultiRaster(cfg_in, enc.S)

   exportgraphics(gcf, ['/home/williamslab/Williams Lab Dropbox/Eric Carmichael/CIHR_2023_September' filesep out.meta.subject '_' out.meta.session  '_raster.pdf'], 'ContentType', 'vector');

%%  REM rad

data_dir = '/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Radial/Radial_ephys/M30/M30_2023-06-06_Rad1'; 

save_dir = '/home/williamslab/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Radial/Inter'; 


if strcmp(computer, 'GLNXA64')
    
    addpath(genpath('/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared'));
    
    addpath(genpath('/home/williamslab/Documents/Github/CEH2'))
    
end

%% 

out = Master_REM_radial_preprocess([], data_dir)

save([save_dir filesep strrep([out.meta.subject '_' out.meta.date '_' out.meta.session], '-', '_') '.mat'],'out');  

%% try it for a session


maps_e = MS_get_place_field([], restrict(out.Encode.S, out.Encode.trials(1,:), out.Encode.trials(2,:)), restrict(out.Encode.pos, out.Encode.trials(1,:), out.Encode.trials(2,:)), out.Encode.speed)


%%
maps_r = MS_get_place_field([], restrict(out.Recall.S, out.Recall.trials(1,:), out.Recall.trials(2,:)), restrict(out.Recall.pos, out.Recall.trials(1,:), out.Recall.trials(2,:)), out.Recall.speed)


%% compare maps across encoding/recall

map_overlap = MS_map_overlap(maps_e, maps_r, 1, {'Encoding', 'Recall'}); 
