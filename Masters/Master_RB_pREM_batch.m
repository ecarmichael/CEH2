%% Batch process RB sleep data


% RB data notes: Sample data for 1 mouse (220).
% Within the main folder are 2 experiment-specific subfolders (each mouse underwent testing under at least 2 different conditions (sometimes 3).
%AC is control data for mouse 220 (no MS GABAergic inhibition at any point of experiment), A is test group (MS GABAergic inhibition selectively during REMs).
% Within each experiment-specific subfolder is the nsc file from the best (largest, most stable recording) and the corresponding hypnogram. 
%(Each datapoint is 1 s, data was scored with 5 s resolution (non-overlapping windows); 
%1 = wake, 2 = NREMs, 3 = REMs, 4 = wake + MS GABAergic inhhibition, 5 = NREMs + MS GABAergic inhibition, 6 = REMs + MS GABAergic inhibition).


%% add codebases

addpath(genpath('/home/williamslab/Documents/Github/CEH2'));
addpath(genpath('/home/williamslab/Documents/Github/vandermeerlab/code-matlab/shared'));

% data_dir = 'J:\Scored_NSC_Data_files\Mouse 220\Experiment 1 (AC)';
% data_dir = 'J:\Scored_NSC_Data_files\Mouse 220\Experiment 2 (A)';
data_dir = '/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/RB_data/RB_Data'; 

% LFP_dir = 'J:\Williams_Lab\Jisoo\LFP data\Jisoo'; 


%% loop over sessions and compute the pREM

cd(data_dir)

sub_list = dir('2*');


for iS = 1:length(sub_list)
    
    cd([data_dir filesep sub_list(iS).name]);
    
    sess_list  = dir;
    
    for iSess = 1:length(sess_list)
        if strcmpi(sess_list(iSess).name, '.') || strcmpi(sess_list(iSess).name, '..')
            continue 
        else
            fprintf('<strong>%s</strong>: extracting pREM ...\n', sess_list(iSess).name)
            cd(sess_list(iSess).name)
%             delete *.png
            RB_get_pREM()
            close all
            cd([data_dir filesep sub_list(iS).name]);
        end
    end
    cd(data_dir)
    
    
end
    
    

