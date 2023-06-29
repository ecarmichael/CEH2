function [CSC, S, pos, evt] = Master_NOL_preprocess(cfg_in, Kilo_dir, OE_dir, DLC_dir, save_dir)
%% MS_preprocess_NOL:  Extracts the position, spike, and CSC data for an OE session. 
%
%
%
%    Inputs: 
%    - cfg_in: [struct]  configuration parameters to overwrite the defaults
%
%    - kilo_dir: [path]  directory with the processed kilosort/phy2
%    output
%
%    -  oE_dir: [path]  directory with the OE continuous files
%
%    - DLC_dir: [path] directory/directories with the DLC outputs. Can be a
%    cell array to process multiple directories for one session. 
%
%    Outputs: 
%    -
%
%
%
%
% EC 2023-06-29   initial version 
%
%
%
%% initialize

cfg_def = [];
cfg_def.conv_fac = [6.6  5.8]; 
cfg = ProcessConfig(cfg_def, cfg_in); 


CSC = [];
S = [];
pos = [];
evt = [];


%%  Get the position data if it is there

if ~isempty(DLC_dir);
    
   [pos, behav] = MS_DLC2TSD(DLC_dir);
   
   % convert pixels to cm
   for ii = 1:length(pos.label); 
       
       if strcmpi(pos.label{ii}, 'HD')
           continue
       end
         pos.data(1,:) = pos.data(1,:)./cfg.conv_fac(1); 
         pos.data(2,:) = pos.data(2,:)./cfg.conv_fac(2);
    
       end
   
   end
   pos.data(3,:) = pos.data(,:)./cfg.conv_fac(1); 
   
      
   behav.position(:,1) = behav.position(:,1)./cfg.conv_fac(1); 
   behav.position(:,2) = behav.position(:,2)./cfg.conv_fac(2); 
   behav.speed(:,1) = behav.speed(:,1)./cfg.conv_fac(1); 

    
    
    
    
    
    
end



