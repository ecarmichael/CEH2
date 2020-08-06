%% JC_sleep_bin_sandbox

time_str = regexp(ms_trk.file_names,'\d*','Match');


ms_trk.time_labels = datestr(datestr([time_str{1} ':' time_str{2} ':' time_str{3}]), 'HH:MM:SS');


%% convert to datenum


trk_time = datevec(ms_trk.time_labels,'HH:MM:SS');
seg_time = datevec(ms_seg.time_labels,'HH:MM:SS');

e_time = etime(seg_time,trk_time);

disp(e_time/60)
%% use values in datestr
bin_s = 20; 
time_bins = -120:20: 120;
time_centers = time_bins(1:end-1)+bin_s/2;
    sessions  = {'LTD1', 'LTD2'}%...};
    for iSess = 1:length(sessions)
        
        
        for iB = 1:length(blocks)
            load(this_block)
            this_time = ms_seg.time2trk;
                
                % find where this_time falls in the time bins. 
                
%                 iT = nearest_idx(this_time,time_centers)
                % 
                if strcmp(blocks{iB}, 'REM')
                    REM_time_vec(iT, end+1) = sum(binary)/length(frames);
                elseif strcmp(blocks{iB}, 'SWS')
                    SWS_time_vec(iT, end+1) = sum(binary)/length(frames);
                end
        end
    end
end




