function All_out = Pipeline_Asmbly(fname,bin_s, fig_dir)
%% Pipeline_Asmbly: provides a wrapper for running assembly and reactivation analyses using calcium data. 

if strcmp(computer, 'GLNXA64')
    
    codebase_dir = '/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared';
    ca_dir = '/home/williamslab/Documents/Github/CEH2';
    oasis_dir = '/home/williamslab/Documents/Github/OASIS_matlab';
    
    code_dir = '/home/williamslab/Documents/Github/Dos-Santos Assembly ICA/Dos-Santos Assembly ICA';
    
    RnR_dir = '/home/williamslab/Documents/Github/RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    data_dir = '/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter';
    
    
    
else
    
    codebase_dir = 'C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared';
    ca_dir = 'C:\Users\ecarm\Documents\GitHub\CEH2';
    oasis_dir = 'C:\Users\ecarm\Documents\GitHub\OASIS_matlab';
    
    code_dir = 'C:\Users\ecarm\Downloads\Dos-Santos Assembly ICA\Dos-Santos Assembly ICA';
    
    RnR_dir = 'C:\Users\ecarm\Documents\GitHub\RnR_methods';
    
    % data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3' %C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\Maze_Ca\inter\pv1254\2021_12_17_pv1254_MZD3';
    data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter';
    
end

restoredefaultpath
c_d = cd; 
cd(oasis_dir)
addpath(genpath(oasis_dir));
oasis_setup

addpath(genpath(ca_dir));
addpath(genpath(codebase_dir))
addpath(genpath(RnR_dir));

addpath(code_dir)

if nargin <2
    save_fig = 0;
elseif nargin == 2
    save_fig = 1;
    if ~exist(fig_dir)
        mkdir(fig_dir)
    end
end
cd(c_d)

rng(123, 'twister')
plot_flag = 1;


%% load the data 

% load('ms_trk.mat')
% load('behav_DLC.mat')
this_sess = fname;

load(this_sess);

%% preprocess beahviour
behav = MS_align_data(behav,ms);

move_idx = behav.speed > 2.5;

%% remove questionable cells

ms_trk = ms;
keep_idx = zeros(1,size(ms_trk.RawTraces,2));
keep_idx(1:floor(size(ms_trk.RawTraces,2)*.66)) = 1;

remove_cell_id = find(~keep_idx);

cfg_rem = [];
cfg_rem.remove_idx = find(~keep_idx);
cfg_rem.data_type = 'RawTraces';
ms_trk_cut = MS_Remove_trace(cfg_rem, ms_trk);

if ~isfield(ms_trk, 'deconv') % get the deconvolved trace if not already present.
    ms_trk_cut = MS_append_deconv(ms_trk_cut, 1);
end
% remove inactive cells

keep_idx = sum(ms_trk_cut.Binary, 1) >0;

remove_cell_id_decon = find(~keep_idx);

cfg_rem = [];
cfg_rem.remove_idx = find(~keep_idx);
cfg_rem.data_type = 'deconv';
ms_trk_cut = MS_Remove_trace(cfg_rem, ms_trk_cut);


% trim the REM data as well
all_binary_post_REM(:, remove_cell_id) = [];
all_binary_post_REM(:, remove_cell_id_decon) = [];

% all_deconv_post_REM 

all_binary_pre_REM(:, remove_cell_id) = [];
all_binary_pre_REM(:, remove_cell_id_decon) = [];

%% load the place information

% dir_parts = strsplit(this_sess, filesep);
% parts = strsplit(dir_parts{end}, '_');
% info.task = parts{2};
% if contains(info.task, 'HATD6')
%     info.task = 'HATDSwitch';
% end
% info.subject = parts{1};


load([info.subject '_' info.session '_PCs.mat'])

place = [];

place.centroids = PCs_properties.peak_loc;

% is it a place cell?
place.is = PCs_properties.isPC;

place.map = PCs_properties.tuning_curve_data';

place.MI = PCs_properties.MI;

place.peak_rate = PCs_properties.peak_rate;



place.centroids(remove_cell_id)= [];
place.centroids(remove_cell_id_decon)= [];

place.is(remove_cell_id)= [];
place.is(remove_cell_id_decon)= [];

place.map(remove_cell_id,:)= [];
place.map(remove_cell_id_decon,:)= [];

place.MI(remove_cell_id)= [];
place.MI(remove_cell_id_decon)= [];

place.peak_rate(remove_cell_id)= [];
place.peak_rate(remove_cell_id_decon)= [];

bin = 3; %80/size(place.map,2);
p_bins = 0:bin:100;
place.p_bins = p_bins(1:end)+bin/2;
% see if there are any anxiety cells

% [~,p_sort] = sort(place.centroids);

%% get the initial assemblies

%     opts.threshold.method = 'binshuffling';
%     opts.Patterns.method = 'ICA';
%     opts.Patterns.number_of_iterations = 500;
%     opts.threshold.permutations_percentile= 95;
%     opts.threshold.number_of_permutations= 500;

A_temp = []; A_prog = []; wake_data = []; wake_tvec = []; 
for ii = length(bin_s):-1:1
    
   [A_temp{ii}, A_proj{ii}, wake_data{ii}, wake_tvec{ii}] = MS_PCA_ICA_only(ms_trk_cut, move_idx, bin_s(ii),[]); 
end


for ii = length(bin_s):-1:1
   fprintf('PCA-ICA detected %.0f assemblies using a %.2fs binsize\n', size(A_temp{ii},2), bin_s(ii)) 
end


%% remove assemblies without positive weights
for iB = length(A_temp):-1:1
    
    [P_temp{iB}, P_proj{iB}, P_pos{iB}] = MS_Asmbly_select(A_temp{iB}, A_proj{iB}, 2); 
    
    fprintf('[%.0f/%.0f = %.0f%%] Assemblies had cells with positive weights (%.2fs binsize)\n',size(P_temp{iB},2),size(A_temp{iB},2),  (size(P_temp{iB},2)/size(A_temp{iB},2))*100, bin_s(iB))
end


%% get the spacial tuning of the assemblies
min_N_place = 3; 

Place_temp = []; Place_proj = []; Place_map = []; 
for iB = length(P_temp):-1:1
    
    [map_out{iB}, place_idx{iB}] = MS_Asmbly_map(P_pos{iB}, place, min_N_place);
    
    Place_map{iB} = map_out{iB}; 
    Place_map{iB}(~place_idx{iB}) = []; 
    
    Place_temp{iB} = P_temp{iB}(:,place_idx{iB}); 
    Place_proj{iB} = P_proj{iB}(:,place_idx{iB}); 

    
  fprintf('[%.0f/%.0f = %.0f%%] Assemblies contained at least %0.0f place cells (%.2fs binsize)\n',size(Place_temp{iB},2),size(A_temp{iB},2),  (size(Place_temp{iB},2)/size(A_temp{iB},2))*100, min_N_place, bin_s(iB))

end

%% get the activation locations on the track;
win_s = 2; 
thresh = 10; 
for iB = length(P_temp):-1:1
    
    [P_loc{iB}] = MS_Asmbly_act_loc(P_proj{iB}, wake_tvec{iB}, behav, win_s, thresh, 2/bin_s(iB));
    
end

%% make a summary plot to visualize the assemblies

for iB = length(P_temp):-1:1
    info.bin = bin_s(iB);
   MS_Asmbly_plot(P_temp{iB}, P_pos{iB}, map_out{iB}, P_loc{iB}, fig_dir, info); 
close all
end

%% Load the REM data and compare it to the wake

for iB = length(bin_s):-1:1
% pre REM

   [REM_pre_proj, REM_pre_data, REM_pre_tvec, REM_pre_R_thresh] = MS_Asmbly_ReAct(all_binary_pre_REM, P_temp{iB} ,ms_trk_cut,  bin_s(iB));



% post REM


end
%% check the correlation and pool 



A_corr = [];
int_tvec = ms.time(1)/1000:min(bin_s):ms.time(end)/1000; 
int_tvec = int_tvec(1:end-1); 

for ii = length(A_prog):-1:1
    
    A_1_int = []; A_2_int = [];
    
    for iA = size(A_prog{ii},1):-1:1
        
        A_1_int = interp1(wake_tvec{ii}, A_prog{ii}(iA,:), int_tvec);
        
        
        
        for jj = length(A_prog):-1:1
            
            for iA2 = size(A_prog{jj},1):-1:1
                
                A_2_int = interp1(wake_tvec{jj}, A_prog{jj}(iA2,:),int_tvec);
                
                
                [this_corr this_p] = corrcoef(A_1_int(~isnan(A_1_int)), A_2_int(~isnan(A_2_int)));
                A_corr(iA,iA2) = this_corr(1,2);
                P_corr(iA, iA2) = this_p(1,2); 
                
            end
            
        end
        
    end
    
end

sum(P_corr(logical(triu(ones(size(P_corr)),1))) < 0.05)/length(P_corr(logical(triu(ones(size(P_corr)),1))))
