%% load some OE data


data_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\CIHR_2025\HF\HF_1_2025-09-05_15-53-34_TFC_REC\Record Node 118';

% 

evts = OE_LoadEvents(); 

csc_list = dir('*.continuous');
csc= []; labels = []; 
for ii = 1:length(csc_list)

    if ii == 1
        [data, tvec, info] = load_open_ephys_data(csc_list(ii).name);

        csc = tsd(tvec, data);
                labels{ii} = info.header.channel; 

    else
        [data, ~, info] = load_open_ephys_data(csc_list(ii).name);
        csc.data =[csc.data, data];  
                labels{ii} = info.header.channel; 

    end
    

end

