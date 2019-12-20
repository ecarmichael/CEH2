function Keys = MS_Write_Keys(data_dir)
%% MS_Write_Exp: this will create an Keys file by taking the information from the data files 
%   and putting it into an exicutable Keys.m script. Yo can then run Keys.m
%   to have a structure with all of the values as a struct
%
%   To Use: put the path to the data folder as the input 'data_dir'.  Or run it without an input 
%   while in the data folder.
%
%   Inputs:
%       - data_dir: [string] path to folder containing experimental data
%           example:
%           MS_Write_Keys('/home/ecarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1')
%
%   Outputs"
%       Note: If no outputs are called then the Keys.m will only be written
%       to the directory but not to the workspace. 
%       - Keys: [strcut]   strcuture containing experimental parameters such as
%           Keys.subject = '537'
%           Keys.date = '2019-12-04'
%           ...
%
% EC - 2019-12-20: initial version based off of STATE_Write_Exp 
%
%
%
%   Keys is based on the ExpKeys format used in the van der Meer lab
%
%
%% initialize

if nargin == 0
    data_dir = cd;
end


%% get the basic information form the name

if isunix
    fname = strsplit(cd, '/');
else
    fname = strsplit(cd, '\');
end


fname = fname{end};
fname = strsplit(fname, '_'); % split into date, time, and ids


subject_id =  fname{3}(1:3);
date_id = fname{1};
time = fname{2}; 
day = str2num(fname{3}(strfind(fname{3}, 'day')+3));
sess_id = fname{3}(end-4:end);

fprintf('\nSession data appears to be Subject %s, on %s, day #%.0f, with task: %s\n', subject_id, date_id, day, sess_id)

% open the file to write
fid = fopen(['M' fname{3} '_Keys.m'], 'w');

%% fill in the consistent information
fprintf(fid, '%% ''This Keys.m was generated using MS_Write_Keysp.m'';\n');

fprintf(fid, ['Keys.version = ' num2str(1) ';\n']);
fprintf(fid, 'Keys.species = ''mouse'';\n');
fprintf(fid, 'Keys.experimenter = ''Eva'';\n');
% fprintf(fid, 'Keys.behavior = ''wheel'';\n');
% fprintf(fid, 'Keys.probe = ''Buz32'';\n');

%% Get the flexible subject information
fprintf(fid, '\n%%Subject information\n');

fprintf(fid, ['Keys.subject = ''' subject_id ''';\n']);
fprintf(fid, ['Keys.date = ''' date_id ''';\n']);
fprintf(fid, ['Keys.session = ''' sess_id ''';\n']);
fprintf(fid, 'Keys.genetics = ''J20'';\n');
fprintf(fid, 'Keys.promoter = ''CAMKII'';\n');
fprintf(fid, 'Keys.reporter = ''GCAMP6f'';\n');
fprintf(fid, 'Keys.target = ''dHC'';\n');
fprintf(fid, 'Keys.age = ''NaN'';\n');
fprintf(fid, 'Keys.task_id = ''linear'';\n');
if strfind(fname{3}, 'base')
    fprintf(fid, 'Keys.task = 0 ;%% did they perform a task 0 = no, 1 = yes\n'); 
else
    fprintf(fid, 'Keys.task = 1 ;%% did they perform a task 0 = no, 1 = yes\n'); 
end
fprintf(fid, 'Keys.task_order = {''baseline'', ''linear'', ''post-baseline''};\n');


fprintf(fid, '\n%%Notes\n');
fprintf(fid, 'Keys.notes = '''';\n');


fprintf(fid, '\n%%Recording details\n');
fprintf(fid, 'Keys.LFP_hemisphere = ''R'';\n');
this_dir = dir;
% get the CSC files here. 
nChan = {}; adChan = [];
for iFile = 1:length(this_dir)
    is_CSC = strfind(this_dir(iFile).name, '.ncs');
    if ~isempty(is_CSC)
        nChan{end+1} = ['CSC' this_dir(iFile).name(is_CSC-1) '.ncs'];
        adChan(end+1) = str2double(this_dir(iFile).name(is_CSC-1));
    end
end

% create the Ad channels
fprintf(fid,'Keys.LFP_ad_channels = {');
for iC = 1:length(adChan)
    if iC ~= length(adChan)
        fprintf(fid, ['''AD_chan' num2str(adChan(iC)) ''',']);
    else    
        fprintf(fid, ['''AD_chan' num2str(adChan(iC)) '''}; %% these are the NLX AD channels that correspond to the CSCs, should be equal but just in case\n']);
    end
end

%same for CSC
fprintf(fid,'Keys.LFP_channels = {');
for iC = 1:length(adChan)
    if iC ~= length(adChan)
        fprintf(fid, ['''' nChan{iC} ''',']);
    else    
        fprintf(fid, ['''' nChan{iC} '''}; %% these are the NLX AD channels that correspond to the CSCs, should be equal but just in case\n']);
    end
end




% fprintf(fid, ['Keys.tetrodeDepths = ' num2str(depth) ';\n']);
fprintf(fid, 'Keys.EMG = ''CSC1.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');


if strcmp(subject_id, '537') 
    fprintf(fid, 'Keys.goodCSC = ''CSC5.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');
    fprintf(fid, 'Keys.goodCSC2 = ''CSC6.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');
    
else
    fprintf(fid, 'Keys.goodCSC = ''CSC5.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');
end
fprintf(fid, 'Keys.quality = ''NaN''; %%0 is poor, 1 means cell faded, 2 means ok, 3 means great!, NaN means not yet filled in\n');

% to determine hemisphere based on mouse and date of switch.

fprintf(fid, 'Keys.Scope_hemisphere = ''L'';\n');
fprintf(fid, 'Keys.Scope_version = ''v3'';\n');




% % laser information
% fprintf(fid, '\n%%Laser information\n');
% fprintf(fid, 'Keys.laser_duration = %.1f; %% in ms.\n', 1);
% fprintf(fid, 'Keys.laser_mW = ''NaN''; %% in mW.\n');
% fprintf(fid, ['Keys.laser_wave = ''' num2str(473) ''';\n']);
% fprintf(fid, 'Keys.laser_type = ''Shanghai Laser'';\n');
% fprintf(fid, 'Keys.laser_shutter = ''yes'';\n');
% fprintf(fid, 'Keys.fibre_cable = %.0d ;%% microns\n', 125);
% fprintf(fid, 'Keys.fibre_probe = %.0d ;%% microns\n', 50);
% fprintf(fid, 'Keys.fibre_NA = %.2f ;%% \n', 0.39);


%% experimental variables: NLX digital I/O
fprintf(fid, '\n%%NLX digital I/O codes\n');
fprintf(fid, 'Keys.event_off = ''TTL Input on AcqSystem1_0 board 0 port 3 value (0x0000).''; %% off signal for all events.  Might need to be updated.\n');


fclose(fid);
disp('Keys written')