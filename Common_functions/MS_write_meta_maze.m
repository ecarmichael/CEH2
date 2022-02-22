function meta = MS_write_meta_maze(data_dir, evt)
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
    evt = [];
elseif nargin <2
    evt = []; 
end

if isempty(evt) && exist('Events.mat','file')
    load('Events.mat'); 
end


%% get the basic information form the name
fname = strsplit(data_dir, filesep);
dir_name = strrep(fname{end},'-','_');

fname = strrep(fname{end},'-','_');
fname = strsplit(fname, '_'); % split into date, time, and ids


subject_id =  fname{4};
date_id = [fname{1} '-' fname{2} '-' fname{3}];
sess_id = fname{end};

fprintf('\n<strong>Session data appears to be Subject %s, on %s, with task: %s</strong>\n', subject_id, date_id,  sess_id)

% open the file to write
fid = fopen([dir_name '_meta.m'], 'w');

%% fill in the consistent information
fprintf(fid, '%% ''This Meta.m was generated using %s.m'';\n', mfilename);

fprintf(fid, ['Meta.version = ' num2str(1) ';\n']);
fprintf(fid, 'Meta.species = ''mouse'';\n');
fprintf(fid, 'Meta.experimenter = ''EC'';\n');
fprintf(fid, 'Meta.behavior = ''%s'';\n', sess_id(1:end-1));
fprintf(fid, 'Meta.probe = ''Miniscope + SE'';\n');

%% Get the flexible subject information
fprintf(fid, '\n%%Subject information\n');

fprintf(fid, ['Meta.subject = ''' subject_id ''';\n']);
fprintf(fid, ['Meta.date = ''' date_id ''';\n']);
fprintf(fid, ['Meta.session = ''' sess_id ''';\n']);
fprintf(fid, 'Meta.genetics = ''%s'';\n', subject_id(1:2));
fprintf(fid, 'Meta.promoter = ''Hsyn'';\n');
fprintf(fid, 'Meta.reporter = ''GCamp6f'';\n');
fprintf(fid, 'Meta.target = ''Ca1'';\n');
fprintf(fid, 'Meta.age = ''NaN'';\n');

fprintf(fid, 'Meta.task_order = {''pre_sleep'', ''W maze'', ''post_sleep''};\n');


fprintf(fid, '\n%%Notes\n');
fprintf(fid, 'Meta.notes = '''';\n');


fprintf(fid, '\n%%Recording details\n');
fprintf(fid, 'Meta.LFP_hemisphere = ''L'';\n');

%same for CSC

% fprintf(fid, ['Meta.tetrodeDepths = ' num2str(depth) ';\n']);
fprintf(fid, 'Meta.EMG = ''CSC1.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');


if strcmp(subject_id, 'pv1254') 
    fprintf(fid, 'Meta.goodCSC = ''CSC6.ncs''; %%this channel was referenced to the skull wire, while all others were locally referenced to optimize spikes\n');
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



%% write the event times for each trial and the type of trial. 


if ~isempty(evt)
    
for ir = length(evt.t{1}):-1:1
    dur(ir) = evt.t{2}(ir) - evt.t{1}(ir);
end
[~, idx] = sort(dur, 'descend'); 
this_rec = idx(3); 

evt_r = restrict(evt, evt.t{1}(this_rec), evt.t{2}(this_rec)); 

    
    % get all the trial starts
    trials.choice = evt_r.t{strcmpi(evt_r.label, 'c')}; 
    trials.start_C = evt_r.t{contains(evt_r.label, 'Box Exited')}((nearest_idx3(trials.choice, evt_r.t{contains(evt_r.label, 'Box Exited')},-1))); 
    
    trials.force_L = evt_r.t{strcmpi(evt_r.label, 'fl')}; 
    trials.start_FL = evt_r.t{contains(evt_r.label, 'Box Exited')}(nearest_idx3(trials.force_L, evt_r.t{contains(evt_r.label, 'Box Exited')},-1)); 
    
    trials.force_R = evt_r.t{strcmpi(evt_r.label, 'fr')}; 
    trials.start_FR = evt_r.t{contains(evt_r.label, 'Box Exited')}((nearest_idx3(trials.force_R, evt_r.t{contains(evt_r.label, 'Box Exited')},-1))); 
    
    trials.rew_L = evt_r.t{strcmpi(evt_r.label, 'L')}; 
    trials.enter_L = evt_r.t{contains(evt_r.label, 'L reward Entered')}((nearest_idx3(trials.rew_L, evt_r.t{contains(evt_r.label, 'L reward Entered')},-1))); 
    
    
    trials.rew_R = evt_r.t{strcmpi(evt_r.label, 'R')}; 
    trials.enter_R = evt_r.t{contains(evt_r.label, 'R reward Entered')}((nearest_idx3(trials.rew_R, evt_r.t{contains(evt_r.label, 'R reward Entered')},-1))); 
    
    
%     trials.box = evt_r.t{strcmpi(evt_r.label, 'b')}; 
%     trials.enter_box = evt_r.t{contains(evt_r.label, 'Box Entered')}((nearest_idx3(trials.box , evt_r.t{contains(evt_r.label, 'Box Entered')}, -1))); 

    trials.box = sort([trials.start_C, trials.start_FL,trials.start_FR]) - 30; 
    
    trials.tstart = sort([trials.start_C, trials.start_FL,trials.start_FR]); 
    trials.tend = sort([trials.rew_L, trials.rew_R]) +5; 
    
    if length(trials.tstart) ~= length(trials.tend)
        if length(trials.tstart)-1 == length(trials.tend)
           warning('tstart (%.0f) tend (%.0f) differ by 1. Assuming incomplete last trial', length(trials.tstart),length(trials.tend)); 
        else
           error('%s: tstart (%.0f) tend (%.0f) differ. Determine why...',mfilename, length(trials.tstart),length(trials.tend)); 
        
        end
    end
    
    % convert to ms times; 
    fnames = fieldnames(trials); 
    for iF = 1:length(fnames)
        trials_ms.(fnames{iF}) = behav.time(nearest_idx3(trials.(fnames{iF}), evt_r.t{end}));
    end
end
% %% test event
% 
% % close(101)
% figure(101)
% subplot(1,2,1)
% pos_r = restrict(pos, evt.t{1}(this_rec), evt.t{2}(this_rec)); 
% hold on
% plot(pos_r.data(1,:), pos_r.data(2,:), '.k')
% 
% plot(pos_r.data(1,nearest_idx3(trials.start_C, pos_r.tvec)), pos_r.data(2,nearest_idx3(trials.start_C, pos_r.tvec)), 'sb', 'markersize', 12, 'markerfacecolor', 'b')
% plot(pos_r.data(1,nearest_idx3(trials.start_FR, pos_r.tvec)), pos_r.data(2,nearest_idx3(trials.start_FR, pos_r.tvec)), 'sr', 'markersize', 12, 'markerfacecolor', 'r')
% plot(pos_r.data(1,nearest_idx3(trials.start_FL, pos_r.tvec)), pos_r.data(2,nearest_idx3(trials.start_FL, pos_r.tvec)), 'sm', 'markersize', 12, 'markerfacecolor', 'm')
% 
% plot(pos_r.data(1,nearest_idx3(trials.enter_R, pos_r.tvec)), pos_r.data(2,nearest_idx3(trials.enter_R, pos_r.tvec)), 'dr', 'markersize', 12, 'markerfacecolor', 'r')
% plot(pos_r.data(1,nearest_idx3(trials.enter_L, pos_r.tvec)), pos_r.data(2,nearest_idx3(trials.enter_L, pos_r.tvec)), 'db', 'markersize', 12, 'markerfacecolor', 'b')
% 
% 
% plot(pos_r.data(1,nearest_idx3(trials.box, pos_r.tvec)), pos_r.data(2,nearest_idx3(trials.box, pos_r.tvec)), 'ob', 'markersize', 12, 'markerfacecolor', 'b')
% 
% % plot(pos_r.data(1,nearest_idx3(evt.t{find(contains(evt.label, 'Box Exited'))}, pos_r.tvec)), pos_r.data(2,nearest_idx3(evt.t{find(contains(evt.label, 'Box Exited'))}, pos_r.tvec)), 'sb', 'markersize', 12, 'markerfacecolor', 'b')
% % plot(pos_r.data(1,nearest_idx3(evt.t{find(contains(evt.label, 'fl'))}, pos_r.tvec)), pos_r.data(2,nearest_idx3(evt.t{find(contains(evt.label, 'fl'))}, pos_r.tvec)), 'sr', 'markersize', 12, 'markerfacecolor', 'r')
% 
% legend({'tracking', 'choise', 'Force R', 'Force L', 'Rew R', 'Rew L', 'Box ent'})
% 
% subplot(1,2,2)
% c_ord = linspecer(length(trials.tend)); 
% hold on
% for ii  = 1:length(trials.tend)
%     this_pos = restrict(pos, trials.tstart(ii), trials.tend(ii));
%     
%     plot(this_pos.data(1,:), this_pos.data(2,:), '.', 'color', c_ord(ii,:)); 
%     pause(0.5)
% end
% 
% %%  event plots using behav. 
% figure(103)
% c_ord = linspecer(length(trials_ms.tend)); 
% hold on
% for ii  = 1:length(trials_ms.tend)
% %     this_pos = nearest_idx3(trials_ms.tstart(ii), behav.time);
%     
%     plot(behav.position(nearest_idx3(trials_ms.tstart(ii), behav.time):nearest_idx3(trials_ms.tend(ii), behav.time),1), behav.position(nearest_idx3(trials_ms.tstart(ii), behav.time):nearest_idx3(trials_ms.tend(ii), behav.time),2), '.', 'color', c_ord(ii,:)); 
%     pause(0.5)
% end
% %%



fclose(fid);
disp('Meta written')