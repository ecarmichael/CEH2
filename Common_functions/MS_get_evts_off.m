function on_off_iv = MS_get_evts_off(evts, pattern)
%% MS_get_evts_off: takes a NLX evts structure and finds the corresponding 'off' for each 'on TTL
%
%
%
%    Inputs: 
%    - evts: [struct]        events structure in the TS format. 
%
%    - patter: [string]      complete or partial label for the 'on' pattern
%
%    Outputs: 
%    - on_off_iv
%
%
%
%
% EC 2023-08-03   initial version 
%
%
%
%% initialize

on_idx = contains(evts.label, pattern);  % find the correct label. 


off_idx = contains(evts.label, [pattern(1:end-9) '(0x0000).']);  % find the corresponding off label


%% toDo checks for patterns that have mulitple off signals. 


% convert to 'iv' (inverval) format for simplicity
on_off_iv = iv(evts.t{on_idx}, evts.t{off_idx}); 
