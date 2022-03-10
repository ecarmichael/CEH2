% sandbox_EV_fig6

% attempt to recreate Fig6 from EV's paper




%% fig6 b 'Assemblies control vs J20+ in NREM and REM

dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\final_analysis\cvsj\1' 


load('assembly_stats.mat');


 rmm = final_stats{2}{end-3}; % get the repeated measures model object. 
 
 % grab the data
NREM = table2array(rmm.BetweenDesign(:,1));
REM= table2array(rmm.BetweenDesign(:,2));
REM_grp = table2array(rmm.BetweenDesign(:,3));
[REM_grp_id, ~,REM_grp_val] = unique(REM_grp); % convert to values; 


control_NREM = NREM(REM_grp_val == 1); 
j20_NREM = NREM(REM_grp_val == 2); 
control_REM = REM(REM_grp_val == 1); 
j20_REM = REM(REM_grp_val == 2); 

figure(101)
bar([mean(control_NREM), mean(j20_NREM);  mean(control_REM), mean(j20_REM)])


%% fig 6 I



 rmm = final_stats{3}{2}; % get the repeated measures model object. 
 
 % grab the data
NREM = table2array(rmm.BetweenDesign(:,1));
REM= table2array(rmm.BetweenDesign(:,2));
REM_grp = table2array(rmm.BetweenDesign(:,3));
[REM_grp_id, ~,REM_grp_val] = unique(REM_grp); % convert to values; 


control_NREM = NREM(REM_grp_val == 1); 
j20_NREM = NREM(REM_grp_val == 2); 
control_REM = REM(REM_grp_val == 1); 
j20_REM = REM(REM_grp_val == 2); 

figure(101)
clf
bar([mean(control_NREM), mean(j20_NREM);  mean(control_REM), mean(j20_REM)])

