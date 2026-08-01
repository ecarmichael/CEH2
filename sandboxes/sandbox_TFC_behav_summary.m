%% sandbox_TFC_behav summary

%% grab all the sessions

s_list = dir('*TFC*'); % get the list of TFC sessions. 

% loop over sessions and collect the arduino outputs. 
for iS = 1:length(s_list)
    
    this_csv = dir(fullfile([s_list(iS).folder filesep s_list(iS).name], '*.csv')); 

    behav_log = readtable([this_csv.folder filesep this_csv.name]); 

    vr = HF_load_VR([this_csv.folder filesep this_csv.name]); 

% get the events files

    evts_fname  = dir(fullfile([this_csv.folder filesep 'Record Node 117'], '*-A.events'));

    evts = OE_LoadEvents([evts_fname.folder filesep evts_fname.name], 30000);


%session info
parts = strsplit(this_csv.folder, filesep);
parts = parts{end};
info.sess = parts(strfind(parts, 'TFC'):end); 
info.sub = parts(1:7); 
info.date = datetime(parts(9:18), 'inputformat', 'yyyy-MM-dd'); 

if strcmpi(info.sess, 'TFCD3')
    HF_TFC_behav_cond([], evts)

end


% end



%% use the evts files instead

