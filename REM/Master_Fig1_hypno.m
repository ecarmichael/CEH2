%% Master_Fig1_hypno_summary


% collects all of the hypnograms and computes basic stats and summary
% plots;

hypno_dir = '/Users/ecar/Williams Lab Dropbox/Eric Carmichael/JisooProject2020/2020_Results_aftercutting/Across_episodes/Inter';

inter_dir = '/Users/ecar/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter/Hypno'; 
if ~exist(inter_dir, "dir"); mkdir(inter_dir); end
%% copy all of the hypnograms (only need to run once)

filelist = dir(fullfile(hypno_dir, ['**' filesep '*Hypno.mat*']));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list

% remove pREM
rm_idx = []; 
for ii = length(filelist):-1:1
    if contains(filelist(ii).folder, 'pREM')
        rm_idx(ii) = 1;
    else
        rm_idx(ii) = 0;
    end
end

filelist(logical(rm_idx)) = []; 

% copy the hypno but named after the session. 

for ii  = 1:length(filelist)
    pv_idx = strfind(lower(filelist(ii).folder), 'pv'); 
    s_name = filelist(ii).folder(pv_idx(end):end); 

    copyfile([filelist(ii).folder filesep filelist(ii).name], [inter_dir filesep 'Hypno_' upper(s_name) '.mat'])

    
end
%%%%%%%%%%%%%%%%
%% Load each hypno and grab all the data


cd(inter_dir)
h_list = dir('Hypno*'); 

for ii = 1:length(h_list)
    load(h_list(ii).name)




end


