function  out = MS_DLC_score_freezing(fname)
%% MS_score_freezing: score freezing based on movement in DLC tracking data. if multiple body parts, average across to get a better measure. 









%% initialize 

% find the label file



disp(9)
pos = MS_DLC2TSD(fname)