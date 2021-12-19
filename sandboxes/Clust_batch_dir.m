%% Run Clust on all sessions

% will check each session in a dir for FD and if one isn't found it will
% run Mxx_Clust_batch on that sessions

main_dir = cd; 
parts = strsplit(cd, filesep);
Sub = parts{end};

dir_list = dir('M26*');

sess_list = {};
for ii = 1:length(dir_list)
    if contains(dir_list(ii).name, Sub)
        sess_list = [sess_list dir_list(ii).name];
    end
end
    

% loop over directories and run corresponding clust_batch script

for iS  = 1:length(sess_list)
    cd(sess_list{iS})

    if exist('FD','dir') == 7 
                cd(main_dir); 

        continue
    elseif  ~isempty(dir('*.ntt'))
    eval([Sub '_Clust_batch']); 
            cd(main_dir); 

    end

end
