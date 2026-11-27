%% sandbox_LFP_pox_preprocess


data_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\PoxR1\LFP'; 


%% loop over sessions

sess_list = dir([data_dir filesep 'Pox*']); 

for ii = 1:length(sess_list)

    cfg_csc.desired_sampling_frequency = 2000; 
    cd([sess_list(ii).folder filesep sess_list(ii).name])
    this_csc = MS_LoadCSC(cfg_csc);


    parts = strsplit(sess_list(ii).name, '_'); 
    sub = parts{1};
    sess = sess_list(ii).name(length(parts{1})+2:end);


    % split based on cohort

    if contains(sub, 'Pox133a')
        %pox133a
        idx = find(contains(this_csc.label, {'CSC2.ncs', 'CSC3.ncs', 'CSC6.ncs', 'CSC7.ncs'})); 
        all_data.Pox133a.(['sess_' sess]).csc = MS_select_tsd(this_csc, idx); 
        all_data.Pox133a.(['sess_' sess]).csc.label = {'Sub1_1', 'Sub2_2', 'CA1_1', 'CA1_2'}; 
        %pox133c
        idx = find(contains(this_csc.label, {'CSC9.ncs', 'CSC10.ncs', 'CSC12.ncs', 'CSC14.ncs' 'CSC15.ncs' 'CSC16.ncs'})); 
        all_data.Pox133c.(['sess_' sess]).csc = MS_select_tsd(this_csc, idx); 
        all_data.Pox133c.(['sess_' sess]).csc.label = {'Sub1_1', 'Sub2_2', 'REF', 'CA1_1', 'CA1_2', 'CA1_3'}; 
        % pox133b
        idx = find(contains(this_csc.label, {'CSC17.ncs', 'CSC18.ncs', 'CSC21.ncs', 'CSC22.ncs', 'CSC23.ncs', 'CSC24.ncs'})); 
        all_data.Pox133b.(['sess_' sess]).csc = MS_select_tsd(this_csc, idx); 
        all_data.Pox133b.(['sess_' sess]).csc.label = {'Sub1_1', 'Sub2_2', 'CA1_1', 'CA1_2', 'CA1_3', 'CA1_4'}; 
    else
        all_data{ii}.cohort = 'KO';
    end
end
