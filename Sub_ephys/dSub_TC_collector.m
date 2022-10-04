%% dSub Lin position summary
clear all; close all
cd('C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\dSubiculum\inter\Spatial_screening_c')
t_files = dir('*tc.mat'); 

TC_R = [];
TC_L = [];
TC_R_z = [];
TC_L_z = [];
for iC = 1:length(t_files)
   load(t_files(iC).name)
    TC = this_TC; 
   for ii = 1:length(TC.TCs)
       if ~isempty(TC.TCs.spk_tc_R_z)
           TC_R_z(end+1,:) = TC.TCs.spk_tc_R_z;
           TC_R(end+1,:) = TC.TCs.spk_tc_R;
       else
           TC_R_z(end+1,:) = nan(1,size(TC_R_z,2));
           TC_R(end+1,:) = nan(1,size(TC_R_z,2));
       end
       if ~isempty(TC.TCs.spk_tc_L_z)
           TC_L_z(end+1,:) = TC.TCs.spk_tc_L_z;
           TC_L(end+1,:) = TC.TCs.spk_tc_L;
       else
           TC_L_z(end+1,:) = nan(1,size(TC_L_z,2));
           TC_L(end+1,:) = nan(1,size(TC_L_z,2));
       end
   end
end

%% sort them based on peak
[~, l_idx] = max(TC_L_z,[],2); 
[~, ls_idx] = sort(l_idx); 
TC_L_z_s = TC_L_z(ls_idx,:); 
TC_L_z_m = TC_L_z_s./max(TC_L_z_s,[],2); 

[~, r_idx] = max(TC_R_z,[],2); 
[~, rs_idx] = sort(l_idx); 
TC_R_z_s = TC_R_z(rs_idx,:); 
TC_R_z_m = TC_R_z_s./max(TC_R_z_s,[],2); 




%% generate a pair of plot

subplot(2,2,1)
title('Left trials')
imagescnan(TC.bin_c,1:size(TC_L_z_s), TC_L_z_s)
set(gca, 'XTick', TC.bin_ticks, 'XTickLabel', TC.bin_tick_labels, 'XTickLabelRotation', 25)
set(gca, 'YTick', 1:35, 'YTickLabel', ls_idx)
vline(TC.bin_ticks_lines(2:end), 'r');
ylabel('Cell ID')

subplot(2,2,2)
title('Right trials')
imagescnan(TC.bin_c,1:size(TC_R_z_s), TC_R_z_s)
set(gca, 'XTick', TC.bin_ticks, 'XTickLabel', TC.bin_tick_labels, 'XTickLabelRotation', 25)
set(gca, 'YTick', 1:35, 'YTickLabel', rs_idx)

vline(TC.bin_ticks_lines(2:end), 'r');


subplot(2,2,3)
title('Left trials')
imagescnan(TC.bin_c,1:size(TC_L_z_m), TC_L_z_m)
set(gca, 'XTick', TC.bin_ticks, 'XTickLabel', TC.bin_tick_labels, 'XTickLabelRotation', 25)
vline(TC.bin_ticks_lines(2:end), 'r');
set(gca, 'YTick', 1:35, 'YTickLabel', ls_idx)
ylabel('Cell ID')

subplot(2,2,4)
title('Right trials')
imagescnan(TC.bin_c,1:size(TC_R_z_m), TC_R_z_m)
set(gca, 'XTick', TC.bin_ticks, 'XTickLabel', TC.bin_tick_labels, 'XTickLabelRotation', 25)
vline(TC.bin_ticks_lines(2:end), 'r');
set(gca, 'YTick', 1:35, 'YTickLabel', rs_idx)

tightfig

set(gcf, 'Position', [381    62   997   910])
