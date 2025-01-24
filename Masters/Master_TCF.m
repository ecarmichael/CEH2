function Master_TCF(data_dir)
%% Master_TCF:  loops over TCF sessions and prepares them for freezing detection based on the protocol. 
%
%
%
%    Inputs: 
%    - data_dir: [path]   directory with the pose files for each session.
%    Assumes they are recoreded as one file per mouse per session. 
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2025-01-22   initial version 
%
%
%% Protocols

% from delpech et al. 2021:
% https://pmc.ncbi.nlm.nih.gov/articles/PMC8763211/#S10. 
% " 4 weeks after the injections in the context “A” and allowed to explore for 240 s,
% at which point a 20-s tone (85 Db, 2000 Hz) was played, followed by a 20-s trace and 
% then a 2-s, 0.75 mA foot shock. This was repeated two more times starting at 402 and 
% 564 s for a total time of 706 s. On day 2 (24 hours later), mice were placed in a 
% context “B” and allowed to explore for 240 s, at which point the same tone as day 1 
% was played for 60 s, followed by 180 s of no-tone (post-tone period). This was repeated 
% two more times for a total time in the context B of 960 s. On day 3 (24 hours later), 
% mice were placed back in context A and allowed to explore for 300 s. Freezing behavior 
% was recorded during all the time spent in either context by an experimenter blind 
% to the genotype and treatment."

TFC1.baseline = [0 240];
TFC1.tone = [240 260; 402 422; 564 584];
TFC1.trace = [260 280; 422 442; 584 604];
TFC1.shock = [280 282; 442 444; 604 606]; 
TFC1.ITI = [282 402; 444 564; 606 706];

TFC2.baseline = [0 240]; 
TFC2.tone = [240 300; 480 540; 720 780];
TFC2.ITI = [300 480; 540 720;780 960]; 

TFC3.baseline = [0 300]; 


%% initialize

f_list = dir([data_dir filesep '*.mp4']); 

rm_idx = zeros(1,length(f_list)); 
for iF = 1:length(f_list)
    if contains(f_list(iF).name, 'labeled.mp4')
        rm_idx(iF) = 1;
    else
        rm_idx(iF) = 0;
    end
end

f_list(logical(rm_idx)) = []; 

for iF = 1:length(f_list)
    fprintf('%s\n',f_list(iF).name)
end


%% load the LED on table

TFC_tab = readtable('CFT_Frame - Sheet1.csv'); 

%% loop over sessons

for iF = 1:length(f_list)
    
    info = [];
    info.subject = f_list(iF).name(strfind(f_list(iF).name, 'Pox'):strfind(f_list(iF).name, '.mp4')-1);
    info.sess = f_list(iF).name(strfind(f_list(iF).name, 'TFC'):strfind(f_list(iF).name, 'TFC')+3);
    info.date = f_list(iF).name(1:strfind(f_list(iF).name, 'TFC')-2);
    
    if strcmp(info.sess, 'TFC1')
        proto = TFC1;
    elseif strcmp(info.sess, 'TFC2')
        proto = TFC2;
    elseif strcmp(info.sess, 'TFC3')
        proto = TFC3;
    end
    
    % get the table info for the lED on frame. 
    this_tab = find(contains(TFC_tab.Subject, info.subject));
    
    out.(info.subject).(info.sess) = MS_DLC_score_freezing(f_list(iF).name,proto, TFC_tab.(info.sess)(this_tab)); 

    
    
    
    
    
    
    
end
    
    
