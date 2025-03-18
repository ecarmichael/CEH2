function [data_out] = CE_starmaze_SWR(data_dir, Chan_to_use)



%% CE_starmaze_swr

cd(data_dir)

csc_dir = dir(fullfile(data_dir, '**', '*.continuous'));

csc_dir = csc_dir(1).folder; 

cd(csc_dir)

d = dir; 
for ii  = 1:length(dir)
    if contains(d(ii).name, 'continuous')
   chan{ii} =str2double(d(ii).name(strfind(d(ii).name, 'CH')+2:strfind(d(ii).name, '.continuous')-1)); 
   rm_idx(ii) = 0; 
    else
        chan{ii} = [];
        rm_idx(ii) = 1; 
    end
end
chan(logical(rm_idx)) = []; 
%%
for ii = 1:22%length(chan)
    cfg = [];
    cfg.fc = {['100_RhythmData-B_CH' num2str(chan{ii}) '.continuous']};
    % cfg.decimateByFactor = 2; % leave empty for reactive decimation using cfg.desired_sampling_frequency
    cfg.desired_sampling_frequency = 2000; % helps with speed.
    
    csc = OE_old_csc2TSD(cfg);
    
    figure(chan{ii})
    plot(csc.tvec, csc.data)
    xlim([ 1038.5 1039])
end
%% get the SWR
% close all
swr = MS_SWR_detector(csc, csc.label{1}, 1)