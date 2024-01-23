function All_out = Pipeline_Asmbly(fname, fig_dir)
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

%% get the initial assemblies
bin_s = [.25, .5]; 

    opts.threshold.method = 'binshuffling';
    opts.Patterns.method = 'ICA';
    opts.Patterns.number_of_iterations = 500;
    opts.threshold.permutations_percentile= 95;
    opts.threshold.number_of_permutations= 500;

A_temp = []; A_prog = []; wake_data = []; wake_tvec = []; 
for ii = length(bin_s):-1:1
    
   [A_temp{ii}, A_prog{ii}, wake_data{ii}, wake_tvec{ii}] = MS_PCA_ICA_only(ms_trk_cut, move_idx, bin_s(ii),[], opts); 
end


for ii = length(bin_s):-1:1
   fprintf('PCA-ICA detected %.0f assemblies using a %.2fs binsize\n', size(A_temp{ii},2), bin_s(ii)) 
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
                
                
                p_corr = corrcoef(A_1_int(~isnan(A_1_int)), A_2_int(~isnan(A_2_int)));
                A_corr(iA,iA2) = p_corr(1,2);
                
            end
            
        end
        
    end
    
end

