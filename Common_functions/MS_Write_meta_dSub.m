function meta = MS_Write_meta_dSub(data_dir)
%% MS_Write_meta_dSub: this will create an Meta file by taking the information from the data files 
%   and putting it into an exicutable Meta.m script. Yo can then run Meta.m
%   to have a structure with all of the values as a struct
%
%   To Use: put the path to the data folder as the input 'data_dir'.  Or run it without an input 
%   while in the data folder.
%
%   Inputs:
%       - data_dir: [string] path to folder containing experimental data
%           example:
%           MS_Write_Meta('/home/ecarmichael/Documents/Williams_Lab/2019-12-04_11-10-01_537day0base1')
%
%   Outputs
%       Note: If no outputs are called then the Meta.m will only be written
%       to the directory but not to the workspace. 
%       - Meta: [strcut]   strcuture containing experimental parameters such as
%           Meta.subject = '537'
%           Meta.date = '2019-12-04'
%           ...
%
% EC - 2019-12-20: initial version based off of STATE_Write_Exp 
%
%
%
%   Meta is based on the ExpMeta format used in the van der Meer lab
%
%
%% initialize

if nargin == 0
    data_dir = cd;
end


%% get the basic information form the name
fname = strsplit(data_dir, filesep);
dir_name = strrep(fname{end},'-','_');

fname = strrep(fname{end},'-','_');
fname = strsplit(fname, '_'); % split into date, time, and ids


subject_id =  fname{1};
date_id = [fname{2} '-' fname{3} '-' fname{4}];
sess_id = fname{end};

fprintf('\n<strong>Session data appears to be Subject %s, on %s, with task: %s</strong>\n', subject_id, date_id,  sess_id)

% open the file to write
fid = fopen([dir_name '_meta.m'], 'w');

%% fill in the consistent information
fprintf(fid, '%% ''This Meta.m was generated using %s.m'';\n', mfilename);

fprintf(fid, ['Meta.version = ' num2str(1) ';\n']);
fprintf(fid, 'Meta.species = ''mouse'';\n');
fprintf(fid, 'Meta.experimenter = ''EC / EP'';\n');
fprintf(fid, 'Meta.behavior = ''OF'';\n');
fprintf(fid, 'Meta.probe = ''Versa 4'';\n');

%% Get the flexible subject information
fprintf(fid, '\n%%Subject information\n');

fprintf(fid, ['Meta.subject = ''' subject_id ''';\n']);
fprintf(fid, ['Meta.date = ''' date_id ''';\n']);
fprintf(fid, ['Meta.session = ''' sess_id ''';\n']);
fprintf(fid, 'Meta.genetics = ''c57blk/6'';\n');
fprintf(fid, 'Meta.promoter = ''N/A'';\n');
fprintf(fid, 'Meta.reporter = ''N/A'';\n');
fprintf(fid, 'Meta.target = ''dSub'';\n');
fprintf(fid, 'Meta.age = ''NaN'';\n');

fprintf(fid, 'Meta.task_order = {''pre_sleep'', ''W maze'', ''OF'', ''post_sleep''};\n');


fprintf(fid, '\n%%Notes\n');
fprintf(fid, 'Meta.notes = '''';\n');


fprintf(fid, '\n%%Recording details\n');
fprintf(fid, 'Meta.LFP_hemisphere = ''R'';\n');

%same for CSC

% fprintf(fid, ['Meta.tetrodeDepths = ' num2str(depth) ';\n']);
fprintf(fid, 'Meta.EMG = ''CSC1.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');


if strcmp(subject_id, 'M21') 
    fprintf(fid, 'Meta.goodCSC = ''CSC8.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');
elseif strcmp(subject_id, 'M23')
    fprintf(fid, 'Meta.goodCSC = ''CSC5.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');
    fprintf(fid, 'Meta.goodCSC2 = ''CSC6.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');
else
    fprintf(fid, 'Meta.goodCSC = ''CSC5.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');
end
fprintf(fid, 'Meta.quality = ''NaN''; %%0 is poor, 1 means cell faded, 2 means ok, 3 means great!, NaN means not yet filled in\n');



%% experimental variables: NLX digital I/O
fprintf(fid, '\n%%NLX digital I/O codes\n');
fprintf(fid, 'Meta.event_off = ''TTL Input on AcqSystem1_0 board 0 port 3 value (0x0000).''; %% off signal for all events.  Might need to be updated.\n');


fclose(fid);
disp('Meta written')