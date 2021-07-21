%% Master_dSub_check
%
% Master script to check cells from the dSub ephys project.  Can be used to
% check any other neuralynx .t data with position files.
%
%

%% initialize some things

addpath(genpath('/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared'))
addpath(genpath('/home/ecarmichael/Documents/GitHub/CEH2'))

inter_dir = '/home/ecarmichael/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/dSubiculum/inter'; % where to save the outputs

data_dir = '/home/ecarmichael/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/dSubiculum/inProcess'; % where to get the data


%% cycle through mice and sessions
sub_list = dir(data_dir);
sub_list(1:2) = [];
for iSub = length(sub_list):-1:1
    subjects{iSub} = sub_list(iSub).name;
end

% loop subjects
for iSub = 1:length(subjects)
    
    % get a list of the good sessions.
    sess_list = dir([data_dir filesep subjects{iSub}]);
    sess_list(1:2) = [];
    sessions = [];
    for iS = length(sess_list):-1:1
        if contains(sess_list(iS).name, 'OF')
            sessions{iS} = sess_list(iS).name;
        end
    end
    sessions(cellfun('isempty', sessions)) = []; 
    
    for iS = 2:-1:1%length(sessions):-1:1
        close all
        % run the screener script saving the output in the inter_dir.
        cd([data_dir filesep subjects{iSub} filesep sessions{iS}])
        dsubscreener_ec(cd, inter_dir);
        
    end
end
