function [REM_stats, SWS_stats] =  MS_compute_REM_prct_JC(Hypno)
%% MS_compute_REM_prct_JC: extract the percentage of REM vs total sleep time in scored data. 
%
%
%
%    Inputs: 
%    - Hypno: [struct]  contains tvec, data (scored values) and labels for
%    corresponding state labels. 

%         Hypno
%         Hypno = 
%           struct with fields:
%             tvec: [35457560×1 double]
%             data: [35457560×1 double]
%             labels: {'Wake'  'SWS'  'REM'  'Quiescence'  'Transition'  'pREM'  'Redo'  'Exit'}
%             cfg: [1×1 struct]
%
%
%
%    Outputs: 
%    - REM_stats: [struct]  contains times and percentages of sleep (REM / (REM + SWS),(REM / length(recording),... ; 
%    - SWS_stats: [struct]  same as REM_stats but for SWS
%
%
%
% EC 2022-01-31   initial version 
%
%
%
%% initialize
SWS_state = find(contains(Hypno.labels,'SWS'));
REM_state = find(contains(Hypno.labels,'REM'));

REM_idx = ismember(Hypno.data, REM_state);
SWS_idx = ismember(Hypno.data, SWS_state);

REM_stats.percent_sleep = (sum(REM_idx) / (sum(REM_idx) + sum(SWS_idx)))*100; 
REM_stats.percent_total = (sum(REM_idx) /length(Hypno.tvec))*100; 
REM_stats.time = sum(REM_idx)*mode(diff(Hypno.tvec)); 

SWS_stats.percent_sleep = (sum(SWS_idx) / (sum(REM_idx) + sum(SWS_idx)))*100; 
SWS_stats.percent_total = (sum(SWS_idx) / length(Hypno.tvec))*100; 
SWS_stats.time = sum(SWS_idx)*mode(diff(Hypno.tvec)); 

fprintf('<strong>%s</strong>: REM <strong>%2.2f%%</strong> of sleep (%.0fsec), <strong>%2.2f%%</strong> of total recording \n',mfilename,  REM_stats.percent_sleep, REM_stats.time, REM_stats.percent_total)
fprintf('<strong>%s</strong>: SWS <strong>%2.2f%%</strong> of sleep (%.0fsec), <strong>%2.2f%%</strong> of total recording\n',mfilename,  SWS_stats.percent_sleep, SWS_stats.time, SWS_stats.percent_total )
