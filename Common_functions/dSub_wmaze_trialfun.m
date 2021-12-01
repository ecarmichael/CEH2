function trials = dSub_wmaze_trialfun(cfg_in, pos); 
%% dSub_wmaze_trialfun: trial function specific to the w_maze. 


%% initialize defaults

cfg_def = [];
cfg_def.type = 'sbox'; % can also be 'continuous'
cfg_def.evt_states = {'R', 'L', 'C', 'S', 'B'}; 
cfg_def.R_Box_xy = [15 40 45 60];
cfg_def.L_Box_xy = [46 40 80 60];
cfg_def.R_Box_ITI_xy = [28 40 36 60];
cfg_def.L_Box_ITI_xy = [55 40 66 60];
cfg_def.R_Arm_xy = [10 20 20 40]; 
cfg_def.L_Arm_xy = [70 20 80 40];
cfg_def.C_Arm_xy = [40 20 55 40];
cfg_def.CL_Inner_Arm_xy = [43 0 60 20];
cfg_def.CR_Inner_Arm_xy = [30 0 43 20];
cfg_def.CL_Outer_Arm_xy = [60 0 80 20];
cfg_def.CR_Outer_Arm_xy = [12 0 30 20];
cfg_def.R_Port_xy = [10 32 20 40];
cfg_def.L_Port_xy = [70 32 80 40];

cfg = ProcessConfig(cfg_def, cfg_in); 

%% get state for each time point. 

states = NaN(size(pos.tvec)); 

segments = {'R_Box_xy', 'L_Box_xy','R_Box_ITI_xy',  'L_Box_ITI_xy', 'R_Arm_xy', 'L_Arm_xy', 'CR_Inner_Arm_xy', 'CL_Inner_Arm_xy',...
    'CR_Outer_Arm_xy', 'CL_Outer_Arm_xy', 'R_Port_xy', 'L_Port_xy','C_Arm_xy'}; 

figure(1)
hold on
scatter(pos.data(1,:), pos.data(2,:), 55, 'k')
c_list = linspecer(length(segments)); 
c_ord(1:2:13,:) = fliplr(c_list(1:7,:));
c_ord(2:2:13,:) = fliplr(c_list(8:end,:));

for iS = 1:length(segments)
seg_idx.(segments{iS}) = (pos.data(1,:) >=cfg.(segments{iS})(1)) & (pos.data(2,:) >=cfg.(segments{iS})(2)) & (pos.data(1,:) <cfg.(segments{iS})(3)) & (pos.data(2,:) <cfg.(segments{iS})(4)); 


states(seg_idx.(segments{iS})) = iS; 

scatter(pos.data(1,seg_idx.(segments{iS})), pos.data(2,seg_idx.(segments{iS})), 55, c_ord(iS, :) )
disp(segments{iS})
pause(.25)


end
% legend(segments)

%% fill in NaNs

states_fill = fillmissing(states, 'previous');

% smooth over 1s window 
states_smooth = smooth(states_fill, round(1/mode(diff(pos.tvec))));

% make it whole to get the sates back to whole numbers.
states_smooth = round(states_smooth);

% plot again

figure(2)
hold on
scatter(pos.data(1,:), pos.data(2,:), 55, 'k')
for iS = 1:length(segments)

    
    scatter(pos.data(1,states_fill == iS), pos.data(2,states_fill == iS), 55, c_ord(iS, :) )
disp(segments{iS})
% pause(1)
end
    axis off
legend(strrep(segments, '_', ' '),'Orientation', 'vertical', 'Location', 'eastoutside')
pause(1);


%% get transitions between states

[trans_val,trans_idx] = findpeaks(abs(diff(states_fill)));%,'MinPeakDistance', round(1/mode(diff(pos.tvec)))); 

states_jitter = states_fill;

for iT = 1:length(trans_idx)
    
%    states_fill(states_fill == find(contains(segments, 'R_Port_xy')))
    
    states_jitter(trans_idx(iT)+1:trans_idx(iT)+30) = states_jitter(trans_idx(iT)+1); 
end







end