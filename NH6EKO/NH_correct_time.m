  % correct the tvec to be in Zeitgeber Time.  Lights on is 08:00 in
     % this experiment with lights off at 20:00. 
        rec_start_time = this_csc{1}.cfg.hdr{1}.TimeCreated(end-5:end); 
        rec_start_time = datetime(datestr([rec_start_time(1:2) ':' rec_start_time(3:4) ':' rec_start_time(5:6)], 'HH:MM:SS'));
        rec_start_time.Format = 'HH:mm:SS'; 
        Z_tvec = all_csc.tvec - all_csc.tvec(1); % correct recording to 0 at start. 
        
        Z_zero = datetime(datestr('07:00:00', 'HH:MM:SS'));
        rec_start = seconds(rec_start_time - Z_zero); 
        Z_tvec = Z_tvec +double(seconds(rec_start_time)); 
        
        Z = seconds(Z_tvec); 
%         Z.Format = 'hh:mm:ss'; 


%%%%%%%
%% sandbox for plots
c_time = all_csc.tvec-all_csc.tvec(1);
plot(c_time, all_csc.data(1,:))
xlim([c_time(1) c_time(end)]); 

h_tick = 0:3600:72*3600; 
set(gca, 'xtick', h_tick, 'xticklabel', 1:72)
