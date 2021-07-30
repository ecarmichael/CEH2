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

% data_dir = %'/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/dSubiculum/inProcess'; % office unix
% inter_dir = %'/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/dSubiculum/inter'; % office unix

% data_dir = 'C:\Users\williamslab\Dropbox (Williams Lab)\Williams Lab Team
% Folder\Eric\dSubiculum\inProcess';
% inter_dir = 'C:\Users\williamslab\Dropbox (Williams Lab)\Williams Lab
% Team Folder\Eric\dSubiculum\inter';

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
    
    for iS = length(sessions):-1:1
        close all
        % run the screener script saving the output in the inter_dir.
        cd([data_dir filesep subjects{iSub} filesep sessions{iS}])
        dsubscreener_ec(cd, inter_dir);
        
    end
end

%% Process all cells

file_list = dir([inter_dir filesep 'All_cells']);

for filename = length(file_list):-1:1
    if contains(file_list(filename).name(1),"M"); 
    names{filename} = file_list(filename).name;
    else
        names{filename}= [];
    end
end
names(cellfun('isempty', names)) = []; 

for iN = length(names):-1:1
    load([inter_dir filesep 'All_cells' filesep names{iN}]);
    pt_ratio(iN)= This_Cell.wave.pt_ratio;
    width(iN)= This_Cell.wave.spike_width;
    fr(iN) = This_Cell.wave.firing_rate;
    This_Cell = [];
end

figure(4000)
scatter3(width, pt_ratio,fr,'MarkerEdgeColor','k','MarkerFaceColor',[0 .75 .75])
xlim([0.1 0.6])
ylim([1 2.5])
xlabel("Width (ms)")
ylabel("Peak/Trough")
zlabel("Firing Rate spikes/s")


%% Run PPC across cells
addpath('/home/ecarmichael/Documents/GitHub/fieldtrip')
ft_defaults;
addpath('/home/ecarmichael/Documents/GitHub/vandermeerlab/code-matlab/shared/ft')

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
    
    for iS = length(sessions):-1:1
        close all
        % run the screener script saving the output in the inter_dir.
        cd([data_dir filesep subjects{iSub} filesep sessions{iS}])
        
        Meta = MS_Load_meta();
        % convert csc to ft format
        data_all = ft_read_neuralynx_interp({Meta.goodCSC});
        
        s_files = dir('*.t'); 
        for iT = length(s_files):-1:1
           spike_ntt = ft_read_spike([s_files(iT).name(1:3) '.ntt']); 
            
           spike = ft_read_spike(s_files(iT).name);
           
           spike.hdr = spike_ntt.hdr; 
           spike.unit{1} = ones(size(spike.timestamp{1}));
           
           data_all = ft_appendspike([], data_all, spike); 
           clear spike spike_ntt
        end
        
        
        MS_get_PPC(cfg_ppc, data_all)
        
    end
end
