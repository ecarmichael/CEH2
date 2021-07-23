
addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));
addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'));

data_dir = 'J:\Williams_Lab\Jisoo\Jisoo_Project\Inter';
Subjects = {'PV1043', 'PV1060','PV1069'};
LFP_dir = 'J:\Williams_Lab\Jisoo\LFP data\Jisoo'; 

for iSub = 1:length(Subjects)
    cd([data_dir filesep Subjects{iSub}])
    sessions = MS_list_dir_names_any(cd, 'T'); % could use MS_list_dir_names(PARAMS.raw_data_dir, {'string'}) to find specific files by replacing 'string' with a thing to find like 'HAT'
    
    for iSess = 1:length(sessions)
        
        cd([data_dir filesep Subjects{iSub} filesep sessions{iSess}])
        if ~exist('ms_resize.mat', 'file')
            continue
        end
%         load('all_binary_pre.mat');
%         tracker{iSub, iSess}.pre = sum(all_binary_pre, 'all');
%         clear all_binary_pre
%         
%         MS_re_binarize_JC(2, cd, cd, 'ms_resize', 'ms_resize');

        save_dir = ['C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\JC_inter\Frequency_mats' filesep sessions{iSess}] ; 
        MS_extract_AMP_Phi(cd, save_dir, 'ms_resize', [], 'estimate');

%         load('all_binary_pre.mat');
%         
%         tracker{iSub, iSess}.post = sum(all_binary_pre, 'all');
%         clear all_binary_pre
    end
end

