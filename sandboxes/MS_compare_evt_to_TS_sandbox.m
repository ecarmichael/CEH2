function evt_out = MS_compare_evt_to_TS_sandbox(cfg_in, evts_in, TS_in)
%% MS_compare_evt_to_TS: used for aligning the Neuralynx .nev events to the Miniscope timestamp.dat 
%   files.  
%
%   WIP: needs ground truth data to validate. 
%
%    Inputs: 
%     - cfg_in: [struct] user configuration paramters. if empty [] then use
%     defaults (see below)
%     - evts_in: [struct] Neuralynx events in the vanderMeer lab format
%     from loadEvents
%     - TS_in: [struct] 
%
%    Outputs: 
%     -
%
%
%
%
% EC 2020-01-14   initial version 
%
%
%% set defaults
cfg_def = [];
cfg_def.t_chan = 3; % which channel in the evt.t to use.
cfg_def.low_val = -10; % where to start on the low end.  Make bigger for more range. 
cfg_def.high_val = 10; % where to start on the high end.  Make bigger for more range. 
cfg_def.check = 0; % use 1 to toggle plotting the output. 
cfg = ProcessConfig2(cfg_def, cfg_in); 

%%
t_idx = 3;

for iRec = 1:length(Rec_idx)
    
    low_val = 10;
    high_val = -10;
    idx_low = NaN; % initialize with something.
    idx_high = NaN;
    
    if iRec == 1
        low_val = 0;  % only move the low_val 'forward' from zero
        while ~isempty(idx_low)
            temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-low_val), evt.t{t_idx}(Rec_idx(iRec+1)));
            idx_low =find(diff(temp_evt.t{t_idx}(1:20)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            %     idx_high =find(diff(temp_evt.t{t_idx}(end-50:end)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            low_val = low_val - 1;
        end
        
        
        while ~isempty(idx_high)
            temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)), evt.t{t_idx}(Rec_idx(iRec+1)-high_val));
            %     idx_low =find(diff(temp_evt.t{t_idx}(1:50)) > mode(diff(temp_evt.t{t_idx}))*1.05)
            idx_high =find(diff(temp_evt.t{t_idx}(end-20:end)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            high_val = high_val +1;
        end
        
        disp(['Corrected indexting ' num2str(0) ' - ' num2str(high_val)])
        
        
        
    elseif iRec == length(Rec_idx)
        high_val = 0; % only move backwards from the end.
        
        while ~isempty(idx_low)
            temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-low_val), evt.t{t_idx}(end));
            idx_low =find(diff(temp_evt.t{t_idx}(1:20)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            %     idx_high =find(diff(temp_evt.t{t_idx}(end-50:end)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            low_val = low_val - 1;
        end
        
        while ~isempty(idx_high)
            temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)), evt.t{t_idx}(end)-high_val);
            %     idx_low =find(diff(temp_evt.t{t_idx}(1:50)) > mode(diff(temp_evt.t{t_idx}))*1.05)
            idx_high =find(diff(temp_evt.t{t_idx}(end-20:end)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            high_val = high_val +1;
        end
        
        
        disp(['Corrected indexting ' num2str(low_val) ' - ' num2str(0)])
        
    else
        
        % loop until you find the index offset that gives no jump in time. for
        % start
        while ~isempty(idx_low)
            temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-low_val), evt.t{t_idx}(Rec_idx(iRec+1)));
            idx_low =find(diff(temp_evt.t{t_idx}(1:15)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            %     idx_high =find(diff(temp_evt.t{t_idx}(end-50:end)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            low_val = low_val - 1;
        end
        
        % loop until you find the index offset that gives no jump in time. for
        % end
        while ~isempty(idx_high)
            temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-low_val), evt.t{t_idx}(Rec_idx(iRec+1)-high_val));
            %     idx_low =find(diff(temp_evt.t{t_idx}(1:50)) > mode(diff(temp_evt.t{t_idx}))*1.05)
            idx_high =find(diff(temp_evt.t{t_idx}(end-15:end)) > mode(diff(temp_evt.t{t_idx}))*1.05);
            high_val = high_val +1;
        end
        
        
        disp(['Corrected indexting ' num2str(low_val) ' - ' num2str(high_val)])
        
    end
    
    times_to_use(iRec,:) = [low_val+1, high_val-1]; % keep the good values for restrictions in the next cell.
    
    if iRec == length(Rec_idx)
        temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-low_val+1), evt.t{t_idx}(end)-high_val-1);
    else
        temp_evt = restrict(evt, evt.t{t_idx}(Rec_idx(iRec)-low_val+1), evt.t{t_idx}(Rec_idx(iRec+1)-high_val-1));
    end
    biggest_jumps(iRec) = max(diff(temp_evt.t{t_idx}));
end

for iRec = 1:length(Rec_idx)
    fprintf('Evt: %.0f Best offsets: start = %.0f  end = %.0f largest jump = %.3f sec\n', iRec, times_to_use(iRec,1), times_to_use(iRec, 2), biggest_jumps(iRec))
    
end