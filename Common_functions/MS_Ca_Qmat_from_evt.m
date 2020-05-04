function Q_mat =  MS_Ca_Qmat_from_evt(cfg_in, ms)
%% MS_Ca_Qmat: generate a Q matrix for the activity (binary) of a calcium channel.  This will bin the Ca activity during specific times.  
%
%
%
%    Inputs: 
%     - cfg_in: [struct]  configurations (see cell below for defaults and
%     options]
%     - ms:     [struct]  miniscope 'ms' common structure.  Should include
%     the Binary field from msExtractBinary_detrendTraces
%     - evts:   [struct]  interval data ('IV') containing the start and
%     stop times for each event.  ie: SWR times.  
%
%    Outputs: 
%     - Q_mat:  [nCell x time] 'TSD' timestamp data with the firing rate in
%     
%
%
%
%
% EC 2020-04-24   initial version 
%
%
%%