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

TFC_tab = readtable('TFC_Frame - Sheet1.csv');

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
        
        out.(info.subject).(info.sess) = MS_DLC_score_freezing(f_list(iF).name,[],proto, TFC_tab.(info.sess)(this_tab), ['figs' filesep info.subject '_' info.sess]);
        
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


%% remove Pox90_M9 for wild behaviour. 

out = rmfield(out, 'Pox90_m9')

%% Resample f_bin 

% t_bin =20;
% 
% 
% for iSub = 1:length(s_list)
%     
%     sess_list = fieldnames(out.(s_list{iSub}));
%     
%     for iSess = 1:length(sess_list)
%         if isempty(out.(s_list{iSub}).(sess_list{iSess}).TFC)
%             continue
%         else
%             
%             f_val = zeros(1, length(0:t_bin:ceil(pos_r.tvec(end))-2));
%             
%             c= 0;
%             for ii = 0:t_bin:ceil(pos_r.tvec(end))-2
%                 c = c+1;
%                 
%                 s_idx = nearest_idx(ii, pos_r.tvec);
%                 
%                 e_idx = nearest_idx(ii+t_bin,pos_r.tvec);
%                 
%                 f_val(c) = sum(fvec(s_idx:e_idx))/length(fvec(s_idx:e_idx));
%                 fprintf('block length %.0f  F: %.2f%%\n', (e_idx - s_idx)*mode(diff(pos_r.tvec)),f_val(c)*100)
%                 
%             end
%             
%         end
%     end


%% sample plot for each subject for both movement and freezing
figure(99)
clf
s_list = fieldnames(out);
c_ord = MS_linspecer(length(s_list)); 

days = {'TFC1', 'TFC2', 'TFC3'}; 

for ii = 1:length(s_list)
    
    for iD = 1:length(days)
        
    subplot(2,3,iD)
    hold on
    plot(out.(s_list{ii}).(days{iD}).pos_r.tvec,out.(s_list{ii}).(days{iD}).pos_r.data(end-1,:), 'color', c_ord(ii,:))

    subplot(2,3,iD+3)
    hold on
    plot(out.(s_list{ii}).(days{iD}).t_bin,out.(s_list{ii}).(days{iD}).f_bin,  'color', c_ord(ii,:))
    end
end


%% collect the outputs

c_ord = MS_linspecer(4); 

s_list = fieldnames(out);
TFC1_out = []; TFC2_out = []; TFC3_out = [];TFC2_bar = [];
TFC2_bar.baseline = [];
TFC2_bar.tone = [];
TFC2_bar.trace = [];
TFC2_bar.ITI = [];
TFC3_bar.baseline = []; 
Geno = []; Subs = []; Sess = []; 
TFC1_geno = []; TFC2_geno = []; TFC3_geno = []; 
TFC2_sub = []; 

for iSub = 1:length(s_list)
    
    sess_list = fieldnames(out.(s_list{iSub}));
    
    for iSess = 1:length(sess_list)
        if isempty(out.(s_list{iSub}).(sess_list{iSess}).TFC)
            continue
        else
                        this_tab = find(contains(TFC_tab.Subject, s_list{iSub}));

            % hold the 60 binned freezing.
            if strcmp(sess_list{iSess}, 'TFC1')
                TFC1_out = [TFC1_out; out.(s_list{iSub}).(sess_list{iSess}).f_bin];
                t_vec = out.(s_list{iSub}).(sess_list{iSess}).t_bin;
                                TFC1_geno(end+1) = TFC_tab.Pox(this_tab); 

            elseif strcmp(sess_list{iSess}, 'TFC2')
                TFC2_out = [TFC2_out; out.(s_list{iSub}).(sess_list{iSess}).f_bin];
                t_vec2 = out.(s_list{iSub}).(sess_list{iSess}).t_bin;
                
                TFC2_bar.baseline(end+1) = mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.baseline);
                TFC2_bar.tone(end+1) =  mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.tone);
                TFC2_bar.trace(end+1) =  mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.trace);
                TFC2_bar.ITI(end+1)=  mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.ITI);
                TFC2_geno(end+1) = TFC_tab.Pox(this_tab);
                TFC2_sub{end+1} = s_list{iSub}; 

            elseif strcmp(sess_list{iSess}, 'TFC3')
                TFC3_out = [TFC3_out; out.(s_list{iSub}).(sess_list{iSess}).f_bin];
                t_vec3 = out.(s_list{iSub}).(sess_list{iSess}).t_bin;
                TFC3_bar.baseline(end+1) = mean(out.(s_list{iSub}).(sess_list{iSess}).TFC.baseline);
                TFC3_geno(end+1) = TFC_tab.Pox(this_tab); 
            end
            Subs{end+1} = s_list{iSub}; 
            Sess{end+1} = sess_list{iSess}; 

            Geno(end+1) = TFC_tab.Pox(this_tab); 
        end
    end
    
end

% simple plots

figure(102)
clf
subplot(2,3,1)
hold on
plot(t_vec/60, TFC1_out)
plot(t_vec/60, mean(TFC1_out), 'k', 'linewidth', 4)
errorbar(t_vec/60, mean(TFC1_out),MS_SEM_vec(TFC1_out)); 
errorbar(t_vec/60, mean(TFC1_out),MS_SEM_vec(TFC1_out)); 

title('Aquisition'); 
ylabel('%freezing')
xlabel('time (s)')

subplot(2,3,2)
hold on
plot(t_vec2/60, TFC2_out, 'linewidth', .5)
plot(t_vec2/60, mean(TFC2_out), 'k', 'linewidth', 4);
errorb(t_vec2/60, mean(TFC2_out), MS_SEM_vec(TFC2_out))

title('Context B + Tone'); 
ylabel('%freezing')
xlabel('time (s)')

subplot(2,3,5)
hold on
b = bar(1:4, [mean(TFC2_bar.baseline), mean(TFC2_bar.tone),mean(TFC2_bar.trace) mean(TFC2_bar.ITI)]);
% b.face
scatter(1+sort(MS_randn_range(length(TFC2_bar.baseline), 1, -.2, .2)), TFC2_bar.baseline,25,  'k', 'filled')
scatter(2+sort(MS_randn_range(length(TFC2_bar.tone), 1, -.2, .2)), TFC2_bar.tone,25,  'k', 'filled')
scatter(3+sort(MS_randn_range(length(TFC2_bar.trace), 1, -.2, .2)), TFC2_bar.trace,25,  'k', 'filled')
scatter(4+sort(MS_randn_range(length(TFC2_bar.ITI), 1, -.2, .2)), TFC2_bar.ITI,25,  'k', 'filled')
xlim([.5 4.5])
set(gca, 'xtick', 1:4, 'XTickLabel', {'Base', 'Tone', 'Trace', 'ITI'})

subplot(2,3,3)
hold on
plot(t_vec3/60, TFC3_out, 'linewidth', .5)
plot(t_vec3/60, mean(TFC3_out), 'k', 'linewidth', 4)


figure(202)
clf
subplot(2,3,1)
hold on
for ii = 1:length(TFC1.tone)
   rectangle('position',  [TFC1.tone(ii,1)/60 0, (TFC1.tone(ii,2) - TFC1.tone(ii,1))/60, 1], 'FaceColor', [.7 .7 .7 .2], 'EdgeColor', 'none')
   rectangle('position',  [TFC1.trace(ii,1)/60 0, (TFC1.trace(ii,2) - TFC1.trace(ii,1))/60, 1], 'FaceColor', c_ord(3,:), 'EdgeColor', 'none')
   rectangle('position',  [TFC1.shock(ii,1)/60 0, (TFC1.shock(ii,2) - TFC1.shock(ii,1))/60, 1], 'FaceColor',[([210 240 10]/255) .8], 'EdgeColor', 'none')

end
plot(t_vec/60, mean(TFC1_out(TFC1_geno == 0,:)), 'color', c_ord(1,:), 'linewidth', 4)
errorbar(t_vec/60, mean(TFC1_out(TFC1_geno == 0,:)),MS_SEM_vec(TFC1_out(TFC1_geno == 0,:)), 'color', c_ord(1,:)); 
errorbar(t_vec/60, mean(TFC1_out(TFC1_geno == 0,:)),MS_SEM_vec(TFC1_out(TFC1_geno == 0,:)), 'color', c_ord(1,:)); 

plot(t_vec/60, mean(TFC1_out(TFC1_geno == 1,:)), 'color', c_ord(2,:), 'linewidth', 4)
errorbar(t_vec/60, mean(TFC1_out(TFC1_geno == 1,:)),MS_SEM_vec(TFC1_out(TFC1_geno == 1,:)), 'color', c_ord(2,:)); 
errorbar(t_vec/60, mean(TFC1_out(TFC1_geno == 1,:)),MS_SEM_vec(TFC1_out(TFC1_geno == 1,:)), 'color', c_ord(2,:)); 

title('Context A: Aquisition'); 
ylabel('%freezing')
xlabel('time (s)')
ylim([0 1])

subplot(2,3,2)
cla
hold on
for ii = 1:length(TFC2.tone)
   rectangle('position',  [TFC2.tone(ii,1)/60 0, (TFC2.tone(ii,2) - TFC2.tone(ii,1))/60, 1], 'FaceColor', [.7 .7 .7 .2], 'EdgeColor', 'none')
   rectangle('position',  [TFC2.trace(ii,1)/60 0, (TFC2.trace(ii,2) - TFC2.trace(ii,1))/60, 1], 'FaceColor', c_ord(3,:), 'EdgeColor', 'none')

end

plot(t_vec2/60, mean(TFC2_out(TFC1_geno == 0,:)), 'color', c_ord(1,:), 'linewidth', 4)
errorbar(t_vec2/60, mean(TFC2_out(TFC1_geno == 0,:)),MS_SEM_vec(TFC2_out(TFC2_geno == 0,:)), 'color', c_ord(1,:)); 
errorbar(t_vec2/60, mean(TFC2_out(TFC1_geno == 0,:)),MS_SEM_vec(TFC2_out(TFC2_geno == 0,:)), 'color', c_ord(1,:)); 

plot(t_vec2/60, mean(TFC2_out(TFC2_geno == 1,:)), 'color', c_ord(2,:), 'linewidth', 4)
errorbar(t_vec2/60, mean(TFC2_out(TFC2_geno == 1,:)),MS_SEM_vec(TFC2_out(TFC2_geno == 1,:)), 'color', c_ord(2,:)); 
errorbar(t_vec2/60, mean(TFC2_out(TFC2_geno == 1,:)),MS_SEM_vec(TFC2_out(TFC2_geno == 1,:)), 'color', c_ord(2,:)); 



title('Context B: Trace Memory'); 
ylabel('%freezing')
xlabel('time (s)')
ylim([0 1])

subplot(2,3,3)
hold on
plot(t_vec3/60, mean(TFC3_out(TFC2_geno == 0,:)), 'color', c_ord(1,:), 'linewidth', 4)
errorbar(t_vec3/60, mean(TFC3_out(TFC2_geno == 0,:)),MS_SEM_vec(TFC3_out(TFC2_geno == 0,:)), 'color', c_ord(1,:)); 
errorbar(t_vec3/60, mean(TFC3_out(TFC2_geno == 0,:)),MS_SEM_vec(TFC3_out(TFC2_geno == 0,:)), 'color', c_ord(1,:)); 

plot(t_vec3/60, mean(TFC3_out(TFC2_geno == 1,:)), 'color', c_ord(2,:), 'linewidth', 4)
errorbar(t_vec3/60, mean(TFC3_out(TFC2_geno == 1,:)),MS_SEM_vec(TFC3_out(TFC2_geno == 1,:)), 'color', c_ord(2,:)); 
errorbar(t_vec3/60, mean(TFC3_out(TFC2_geno == 1,:)),MS_SEM_vec(TFC3_out(TFC2_geno == 1,:)), 'color', c_ord(2,:)); 

title('Context A'': Context Memory'); 
ylabel('%freezing')
xlabel('time (s)')
ylim([0 1])

% plot(t_vec3/60, TFC3_out, 'linewidth', .5)
% plot(t_vec3/60, mean(TFC3_out), 'k', 'linewidth', 4)

subplot(2,3,4)
cla
hold on
% 
% scatter(1+sort(MS_randn_range(length(TFC2_bar.tone(TFC2_geno == 0)), 1, -.1, .1)), TFC2_bar.tone(TFC2_geno == 0),25,  c_ord(1,:), 'filled')
% scatter(2+sort(MS_randn_range(length(TFC2_bar.tone(TFC2_geno == 1)), 1, -.1, .1)), TFC2_bar.tone(TFC2_geno == 1),25,  c_ord(2,:), 'filled')
% 
% scatter(5+sort(MS_randn_range(length(TFC2_bar.trace(TFC2_geno == 0)), 1, -.1, .1)), TFC2_bar.trace(TFC2_geno == 0),25,  c_ord(1,:), 'filled')
% scatter(6+sort(MS_randn_range(length(TFC2_bar.trace(TFC2_geno == 1)), 1, -.1, .1)), TFC2_bar.trace(TFC2_geno == 1),25,  c_ord(2,:), 'filled')

[hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.baseline(TFC2_geno == 0), TFC2_bar.baseline(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [1 2]);

hb(1).FaceColor = 'none';
hb(1).EdgeColor = 'k';

[hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.tone(TFC2_geno == 0), TFC2_bar.tone(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [4 5]);
hb(1).FaceColor = 'none';
hb(1).EdgeColor = 'k';

[hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.trace(TFC2_geno == 0), TFC2_bar.trace(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [7 8]);
hb(1).FaceColor = 'none';
hb(1).EdgeColor = 'k';

[hb, eb, sc, p, stats] = MS_bar_w_err(TFC2_bar.ITI(TFC2_geno == 0), TFC2_bar.ITI(TFC2_geno == 1), [c_ord(1,:);c_ord(2,:)],1,  'ttest2', [10 11]);
hb(1).FaceColor = 'none';
hb(1).EdgeColor = 'k';


y_l = ylim; 
ylim([y_l(1) 1.1])

% text(1.5, y_l(2)*1.1, 'Tone', 'HorizontalAlignment','center','VerticalAlignment','top')

% text(5.5,y_l(2)*1.1, 'Trace', 'HorizontalAlignment','center','VerticalAlignment','top')

% set(gca, 'xtick', [1:2 5:6], 'xticklabels', {'Tau -', 'Tau +', 'Tau -', 'Tau +'}, 'XTickLabelRotation', 45)
set(gca, 'xtick', 1.5:3:11.5, 'XTickLabel', {'Base', 'Tone', 'Trace', 'ITI'})

% y_l = ylim; 



% title('Context B')
ylabel('Trace freezing (%)')


subplot(2,3,5)
hold on
b = bar(1:4, [mean(TFC2_bar.baseline(TFC2_geno == 0)), mean(TFC2_bar.tone(TFC2_geno == 0)),mean(TFC2_bar.trace(TFC2_geno == 0)) mean(TFC2_bar.ITI(TFC2_geno == 0));...
mean(TFC2_bar.baseline(TFC2_geno == 1)), mean(TFC2_bar.tone(TFC2_geno == 1)),mean(TFC2_bar.trace(TFC2_geno == 1)) mean(TFC2_bar.ITI(TFC2_geno == 1))]);
b(1).FaceColor =  [1 1 1]; 
b(2).FaceColor =   [1 1 1]; c_ord(2,:);

scatter(.85+sort(MS_randn_range(length(TFC2_bar.baseline(TFC2_geno == 0)), 1, -.05, .05)), TFC2_bar.baseline(TFC2_geno == 0),25,   c_ord(1,:), 'filled')
scatter(1.85+sort(MS_randn_range(length(TFC2_bar.tone(TFC2_geno == 0)), 1, -.05, .05)), TFC2_bar.tone(TFC2_geno == 0),25,   c_ord(1,:), 'filled')
scatter(2.85+sort(MS_randn_range(length(TFC2_bar.trace(TFC2_geno == 0)), 1, -.05, .05)), TFC2_bar.trace(TFC2_geno == 0),25,   c_ord(1,:), 'filled')
scatter(3.85+sort(MS_randn_range(length(TFC2_bar.ITI(TFC2_geno == 0)), 1, -.05, .05)), TFC2_bar.ITI(TFC2_geno == 0),25,   c_ord(1,:), 'filled')

scatter(1.15+sort(MS_randn_range(length(TFC2_bar.baseline(TFC2_geno == 1)), 1, -.05, .05)), TFC2_bar.baseline(TFC2_geno == 1),25,   c_ord(2,:), 'filled')
scatter(2.15+sort(MS_randn_range(length(TFC2_bar.tone(TFC2_geno == 1)), 1, -.05, .05)), TFC2_bar.tone(TFC2_geno == 1),25,  c_ord(2,:), 'filled')
scatter(3.15+sort(MS_randn_range(length(TFC2_bar.trace(TFC2_geno == 1)), 1, -.05, .05)), TFC2_bar.trace(TFC2_geno == 1),25,  c_ord(2,:), 'filled')
scatter(4.15+sort(MS_randn_range(length(TFC2_bar.ITI(TFC2_geno == 1)), 1, -.05, .05)), TFC2_bar.ITI(TFC2_geno == 1),25,   c_ord(2,:), 'filled')
xlim([.5 4.5])
set(gca, 'xtick', 1:4, 'XTickLabel', {'Base', 'Tone', 'Trace', 'ITI'})
ylim(y_l)





subplot(2,3,6)
hold on
scatter(1+sort(MS_randn_range(length(TFC3_bar.baseline(TFC3_geno == 0)), 1, -.1, .1)), TFC3_bar.baseline(TFC3_geno == 0),25,  c_ord(1,:), 'filled')
scatter(2+sort(MS_randn_range(length(TFC3_bar.baseline(TFC3_geno == 1)), 1, -.1, .1)), TFC3_bar.baseline(TFC3_geno == 1),25,  c_ord(2,:), 'filled')

[hb, h, p] = MS_bar_w_err(TFC3_bar.baseline(TFC3_geno == 0), TFC3_bar.baseline(TFC3_geno == 1), [.2 .2 .2],1,  'ttest2');

set(gca, 'xtick', 1:2, 'xticklabels', {'Tau -', 'Tau +'}, 'XTickLabelRotation', 45)
ylim([y_l(1) y_l(2)*1.1])
hb(1).FaceColor = 'none';
hb(1).EdgeColor = 'k';
xlim([0 6])

ylim(y_l)
% title('Context A Re-exposure')
ylabel('Contextual freezing (%)')


%% repeated measures
TFC1_vec = []; 
TFC2_vec = []; 

for iSub = length(s_list):-1:1
    
    TFC1_vec(iSub,:) = out.(s_list{iSub}).TFC1.f_bin; 

    
    TFC2_vec(iSub,:) = out.(s_list{iSub}).TFC2.f_bin; 
    
end

fprintf('\n')
for ii = 1:size(TFC1_vec,2)
    fprintf('''t%d'', ', ii);
end
fprintf('\n')

% TFC1

tfc1_tab = array2table([TFC1_geno', TFC1_vec], 'VariableNames', {'Geno',...
't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't9', 't10', 't11', 't12', 't13', 't14', 't15', 't16', 't17', 't18', 't19', 't20', 't21', 't22', 't23', 't24',...
});

tfc1_tab.Geno = categorical(tfc1_tab.Geno); 

rm = fitrm(tfc1_tab, 't1-t24 ~ Geno', 'WithinDesign',  out.(s_list{1}).TFC1.t_bin);

display('TFC1')
ranova_tbl = ranova(rm)



% TFC2

tfc2_tab = array2table([TFC2_geno', TFC2_vec], 'VariableNames', {'Geno',...
    't1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't9', 't10', 't11', 't12', 't13', 't14', 't15', 't16', 't17', 't18', 't19', 't20', 't21', 't22', 't23', 't24', 't25', 't26', 't27', 't28', 't29', 't30', 't31', 't32'...
    });
tfc2_tab.Geno = categorical(tfc2_tab.Geno); 

rm = fitrm(tfc2_tab, 't1-t32 ~ Geno', 'WithinDesign',  out.(s_list{1}).TFC2.t_bin);
display('TFC2')
ranova_tbl = ranova(rm)

%% convert data into csv

table_out = table(TFC2_sub', TFC2_geno',TFC2_bar.baseline',TFC2_bar.tone', TFC2_bar.trace', TFC2_bar.ITI', TFC3_bar.baseline',...
    'VariableNames', {'Subject', 'Genotype', 'TFC2_baseline', 'TFC2_tone', 'TFC2_trace', 'TFC2_ITI','TFC3_baseline'});

writetable(table_out, 'TFC_bar_data.csv')

%% save

set(gcf,'PaperPositionMode','auto');         
set(gcf,'PaperOrientation','landscape');
set(gcf,'Position',[50 50 1200 800]);
print(gcf, '-dpdf',[cd filesep 'figs' filesep 'TFC_summary.pdf'])

end
