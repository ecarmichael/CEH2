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
TFC2.trace = [300 320; 540 560; 780 800] ; 
TFC2.ITI = [320 480; 560 720; 800 960];

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
out = [];

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
    if isempty(this_tab)
        continue
    end
    
    if ~isempty(this_tab)%+ ~isnan(TFC_tab.(info.sess)(this_tab))) == 0
        
        out.(info.subject).(info.sess) = MS_DLC_score_freezing(f_list(iF).name,[],proto, TFC_tab.(info.sess)(this_tab), ['Figs' filesep info.subject '_' info.sess]);
        
    else
        out.(info.subject).(info.sess).out = [];
        out.(info.subject).(info.sess).out.fvec = [];
        out.(info.subject).(info.sess).out.f_bin = [];
        out.(info.subject).(info.sess).out.t_bin = [];
        out.(info.subject).(info.sess).out.TFC = [];
    end
    
    % hold the 60 binned freezing.
    %     if strcmp(info.sess, 'TFC1')
    %         TFC1_out = [TFC1_out, out.(info.subject).(info.sess).f_bin];
    %     elseif strcmp(info.sess, 'TFC2')
    %         TFC2_out = [TFC2_out, out.(info.subject).(info.sess).f_bin];
    %     elseif strcmp(info.sess, 'TFC3')
    %         TFC3_out = [TFC3_out, out.(info.subject).(info.sess).f_bin];
    %     end
    %
end
%




%% Resample f_bin 

t_bin =20;


for iSub = 1:length(s_list)
    
    sess_list = fieldnames(out.(s_list{iSub}));
    
    for iSess = 1:length(sess_list)
        if isempty(out.(s_list{iSub}).(sess_list{iSess}).TFC)
            continue
        else
            
            f_val = zeros(1, length(0:t_bin:ceil(pos_r.tvec(end))-2));
            
            c= 0;
            for ii = 0:t_bin:ceil(pos_r.tvec(end))-2
                c = c+1;
                
                s_idx = nearest_idx(ii, pos_r.tvec);
                
                e_idx = nearest_idx(ii+t_bin,pos_r.tvec);
                
                f_val(c) = sum(fvec(s_idx:e_idx))/length(fvec(s_idx:e_idx));
                %fprintf('block length %.0f  F: %.2f%%\n', (e_idx - s_idx)*mode(diff(pos_r.tvec)),f_val(c)*100)
                
            end
            
        end
    end


%% collect the outputs

s_list = fieldnames(out);
TFC1_out = []; TFC2_out = []; TFC3_out = [];TFC2_bar = [];
TFC2_bar.baseline = [];
TFC2_bar.tone = [];
TFC2_bar.trace = [];
TFC2_bar.ITI = [];
TFC3_bar.baseline = []; 

for iSub = 1:length(s_list)
    
    sess_list = fieldnames(out.(s_list{iSub}));
    
    for iSess = 1:length(sess_list)
        if isempty(out.(s_list{iSub}).(sess_list{iSess}).TFC)
            continue
        else
            % hold the 60 binned freezing.
            if strcmp(sess_list{iSess}, 'TFC1')
                TFC1_out = [TFC1_out; out.(s_list{iSub}).(sess_list{iSess}).f_bin];
                t_vec = out.(s_list{iSub}).(sess_list{iSess}).t_bin;
            elseif strcmp(sess_list{iSess}, 'TFC2')
                TFC2_out = [TFC2_out; out.(s_list{iSub}).(sess_list{iSess}).f_bin];
                t_vec2 = out.(s_list{iSub}).(sess_list{iSess}).t_bin;
                
                TFC2_bar.baseline(end+1) = mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.baseline);
                TFC2_bar.tone(end+1) =  mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.tone);
                TFC2_bar.trace(end+1) =  mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.trace);
                TFC2_bar.ITI(end+1)=  mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.ITI);
                
            elseif strcmp(sess_list{iSess}, 'TFC3')
                TFC3_out = [TFC3_out; out.(s_list{iSub}).(sess_list{iSess}).f_bin];
                t_vec3 = out.(s_list{iSub}).(sess_list{iSess}).t_bin;
                TFC3_bar.baseline(end+1) = mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.baseline);

            end
            
        end
    end
    
end

% simple plots

figure(102)
clf
subplot(2,3,1)
hold on
plot(t_vec/60, mean(TFC1_out), 'k', 'linewidth', 4)
errorb(t_vec/60, TFC1_out, MS_SEM_vec(TFC1_out))

title('Aquisition'); 
ylabel('%freezing')
xlabel('time (s)')

subplot(2,3,2)
hold on
plot(t_vec2/60, TFC2_out, 'linewidth', .5)
plot(t_vec2/60, mean(TFC2_out), 'k', 'linewidth', 4)
title('Context B + Tone'); 
ylabel('%freezing')
xlabel('time (s)')

subplot(2,3,5)
hold on
bar(1:4, [mean(TFC2_bar.baseline), mean(TFC2_bar.tone),mean(TFC2_bar.trace) mean(TFC2_bar.ITI)])
scatter(1+sort(MS_randn_range(length(TFC2_bar.baseline), 1, -.2, .2)), TFC2_bar.baseline,25,  'b', 'filled')
scatter(2+sort(MS_randn_range(length(TFC2_bar.tone), 1, -.2, .2)), TFC2_bar.tone,25,  'b', 'filled')
scatter(3+sort(MS_randn_range(length(TFC2_bar.trace), 1, -.2, .2)), TFC2_bar.trace,25,  'b', 'filled')
scatter(4+sort(MS_randn_range(length(TFC2_bar.ITI), 1, -.2, .2)), TFC2_bar.ITI,25,  'b', 'filled')
xlim([.5 4.5])
set(gca, 'xtick', 1:4, 'XTickLabel', {'Base', 'Tone', 'Trace', 'ITI'})

subplot(2,3,3)
hold on
plot(t_vec3/60, TFC3_out, 'linewidth', .5)
plot(t_vec3/60, mean(TFC3_out), 'k', 'linewidth', 4)

end
