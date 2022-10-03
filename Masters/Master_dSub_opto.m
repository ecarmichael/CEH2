%% sandbox dSub opto stats
Data = readtable('dSub_CMK2_opto - Sheet1.csv');

A0_nt = str2double(Data{4,2:end});
A0_c = str2double(Data{5,2:end});
A0_dur = str2double(Data{7,2:end});

Ab_nt = str2double(Data{10,2:end});
Ab_c = str2double(Data{11,2:end});
Ab_dur = str2double(Data{13,2:end});

B0_nt = str2double(Data{16,2:end});
B0_c = str2double(Data{17,2:end});
B0_dur = str2double(Data{19,2:end});

Bb_nt = str2double(Data{22,2:end});
Bb_c = str2double(Data{23,2:end});
Bb_dur = str2double(Data{25,2:end});

C0_nt = str2double(Data{28,2:end});
C0_c = str2double(Data{29,2:end});
C0_dur = str2double(Data{31,2:end});

Cb_nt = str2double(Data{34,2:end});
Cb_c = str2double(Data{35,2:end});
Cb_dur = str2double(Data{37,2:end});

D0_nt = str2double(Data{40,2:end});
D0_c = str2double(Data{41,2:end});
D0_dur = str2double(Data{43,2:end});

VTEs = readtable('dSub_CMK2_opto - VTEs.csv');

day = str2double(VTEs{2,2:end});
type = VTEs{3,2:end};
A0_vte = str2double(VTEs{4,2:end});
Ab_vte = str2double(VTEs{5,2:end});
B0_vte = str2double(VTEs{6,2:end});
Bb_vte = str2double(VTEs{7,2:end});
C0_vte = str2double(VTEs{8,2:end});
Cb_vte = str2double(VTEs{9,2:end});
D0_vte = str2double(VTEs{10,2:end});

%% get means and SEMs
ArchT_nt = [B0_nt ; Bb_nt ; C0_nt ; Cb_nt];
ArchT_c = [B0_c./B0_nt ; Bb_c./Bb_nt ; C0_c./C0_nt ; Cb_c./Cb_nt]*100;
ArchT_dur = [B0_dur ; Bb_dur ; C0_dur ; Cb_dur];

ctrl_nt = [A0_nt; Ab_nt; D0_nt];
ctrl_c = [A0_c./A0_nt; Ab_c./Ab_nt; D0_c./D0_nt]*100;
ctrl_dur = [A0_dur; Ab_dur; D0_dur];

ArchT_nt_u = mean(ArchT_nt); 
ctrl_nt_u = mean(ctrl_nt); 

ArchT_nt_sem = std(ArchT_nt)./sqrt(length(ArchT_nt)); 
ctrl_nt_sem = std(ctrl_nt)./sqrt(length(ctrl_nt)); 

ArchT_c_u = mean(ArchT_c); 
ctrl_c_u = mean(ctrl_c); 

ArchT_c_sem = std(ArchT_c)./sqrt(length(ArchT_c)); 
ctrl_c_sem = std(ctrl_c)./sqrt(length(ctrl_c)); 

%% VTE get means and SEMs

ArchT_vte = [B0_vte ; Bb_vte; C0_vte ; Cb_vte];
ctrl_vte = [A0_vte ; Ab_vte; D0_vte]; 

ArchT_vte_u = mean(ArchT_vte); 
ctrl_vte_u = mean(ctrl_vte); 

ArchT_vte_sem = std(ArchT_vte)./sqrt(length(ArchT_vte)); 
ctrl_vte_sem = std(ctrl_vte)./sqrt(length(ctrl_vte)); 

% VTE rate

ArchT_vte_r = [B0_vte./B0_nt ; Bb_vte./Bb_nt; C0_vte./C0_nt ; Cb_vte./Cb_nt]*100;
ctrl_vte_r = [A0_vte./A0_nt ; Ab_vte./Ab_nt; D0_vte./D0_nt]*100; 

ArchT_vte_r_u = mean(ArchT_vte_r); 
ctrl_vte_r_u = mean(ctrl_vte_r); 

ArchT_vte_r_sem = std(ArchT_vte_r)./sqrt(length(ArchT_vte_r)); 
ctrl_vte_r_sem = std(ctrl_vte_r)./sqrt(length(ctrl_vte_r)); 
%%  error bar plot for VTE
green = [60 200 50]/255;
c_ord = linspecer(12); 
figure(1)
title('VTE per session')
cla
hold on
rectangle('Position', [16, floor(min([ArchT_vte_u ctrl_vte_u])), 11, ceil(max([ArchT_vte_u ctrl_vte_u])- min([ArchT_vte_u ctrl_vte_u])+1)], 'FaceColor', [.6 .6 .6 .2], 'EdgeColor',[.6 .6 .6 .2])
rectangle('Position', [29, floor(min([ArchT_vte_u ctrl_vte_u])), 7, ceil(max([ArchT_vte_u ctrl_vte_u])- min([ArchT_vte_u ctrl_vte_u])+1)], 'FaceColor', [c_ord(7,:) .2], 'EdgeColor',[c_ord(7,:) .2])

h = errorbar(day, ArchT_vte_u, ArchT_vte_sem, 'color', green, 'LineWidth', 2);
h2 = plot(day, ArchT_vte_u, 's',  'MarkerSize', 14, 'MarkerFaceColor',  green, 'MarkerEdgeColor',  green);

hc = errorbar(day, ctrl_vte_u, ctrl_vte_sem, 'color', [.7 .7 .7], 'LineWidth', 2);
hc2 = plot(day, ctrl_vte_u, 'o',  'MarkerSize', 14, 'MarkerFaceColor',  [.7 .7 .7], 'MarkerEdgeColor',  [.7 .7 .7]);

ylabel('number of trials with VTEs')
xlabel('session')
xlim([0 35])
xline(15,'--k')
xline(28,'--k')
text(0,floor(min([ArchT_vte_u ctrl_vte_u])),  'Transparent Box ''on''', 'VerticalAlignment', 'bottom', 'FontSize', 12)
text(16,floor(min([ArchT_vte_u ctrl_vte_u])),  'Dark Box ''on''', 'VerticalAlignment', 'bottom', 'FontSize', 12)
text(29,floor(min([ArchT_vte_u ctrl_vte_u])),  'Start Arm ''on''', 'VerticalAlignment', 'bottom',  'FontSize', 12)
legend([h2, hc2], {'ArchT(n=4)', 'Control (n=3)'}, 'FontSize', 16, 'Location', 'northwest', 'Box', 'off')

SetFigure([], gcf)

%%  error bar plot for VTE rate
green = [60 200 50]/255;
c_ord = linspecer(12); 
figure(2)
% title('VTE per session')
cla
hold on
rectangle('Position', [16, 0, 11, 100], 'FaceColor', [.6 .6 .6 .2], 'EdgeColor',[.6 .6 .6 .2])
rectangle('Position', [29, 0, 7, 100], 'FaceColor', [c_ord(7,:) .2], 'EdgeColor',[c_ord(7,:) .2])

h = errorbar(day, ArchT_vte_r_u, ArchT_vte_r_sem, 'color', green, 'LineWidth', 2);
h2 = plot(day, ArchT_vte_r_u, 's',  'MarkerSize', 14, 'MarkerFaceColor',  green, 'MarkerEdgeColor',  green);

hc = errorbar(day, ctrl_vte_r_u, ctrl_vte_r_sem, 'color', [.7 .7 .7], 'LineWidth', 2);
hc2 = plot(day, ctrl_vte_r_u, 'o',  'MarkerSize', 14, 'MarkerFaceColor',  [.7 .7 .7], 'MarkerEdgeColor',  [.7 .7 .7]);

ylabel('% of trials with VTEs')
xlabel('session')
xlim([0 35])
xline(15,'--k')
xline(28,'--k')
legend([h2, hc2], {'ArchT(n=4)', 'Control (n=3)'}, 'FontSize', 16, 'Location', 'northeast', 'Box', 'off')

text(0,0,  'Transparent Box ''on''', 'VerticalAlignment', 'bottom', 'FontSize', 18)
text(16,0,  'Dark Box ''on''', 'VerticalAlignment', 'bottom', 'FontSize', 18)
text(29,0,  'Start Arm ''on''', 'VerticalAlignment', 'bottom',  'FontSize', 18)
ylim([0 100])

SetFigure([], gcf)

figure(3)
cla
hold on

h = errorbar(day(1:14), ArchT_vte_r_u(1:14), ArchT_vte_r_sem(1:14), 'color', green, 'LineWidth', 2);
h2 = plot(day(1:14), ArchT_vte_r_u(1:14), 's',  'MarkerSize', 14, 'MarkerFaceColor',  green, 'MarkerEdgeColor',  green);

hc = errorbar(day(1:14), ctrl_vte_r_u(1:14), ctrl_vte_r_sem(1:14), 'color', [.7 .7 .7], 'LineWidth', 2);
hc2 = plot(day(1:14), ctrl_vte_r_u(1:14), 'o',  'MarkerSize', 14, 'MarkerFaceColor',  [.7 .7 .7], 'MarkerEdgeColor',  [.7 .7 .7]);

ylabel('% of trials with VTEs')
xlabel('session')
xlim([0 14.5])
legend([h2, hc2], {'ArchT(n=4)', 'Control (n=3)'}, 'FontSize', 16, 'Location', 'northeast', 'Box', 'off')
ylim([0 100])
text(0,0,  'Transparent Box ''on''', 'VerticalAlignment', 'bottom', 'FontSize', 18)

SetFigure([], gcf)

figure(4)
cla
hold on
rectangle('Position', [14, 0, 14, 100], 'FaceColor', [.6 .6 .6 .2], 'EdgeColor',[.6 .6 .6 .2])

h = errorbar(day(15:27), ArchT_vte_r_u(15:27), ArchT_vte_r_sem(15:27), 'color', green, 'LineWidth', 2);
h2 = plot(day(15:27), ArchT_vte_r_u(15:27), 's',  'MarkerSize', 14, 'MarkerFaceColor',  green, 'MarkerEdgeColor',  green);

hc = errorbar(day(15:27), ctrl_vte_r_u(15:27), ctrl_vte_r_sem(15:27), 'color', [.7 .7 .7], 'LineWidth', 2);
hc2 = plot(day(15:27), ctrl_vte_r_u(15:27), 'o',  'MarkerSize', 14, 'MarkerFaceColor',  [.7 .7 .7], 'MarkerEdgeColor',  [.7 .7 .7]);

ylabel('% of trials with VTEs')
xlabel('session')
xlim([15.5 27.5])
legend([h2, hc2], {'ArchT(n=4)', 'Control (n=3)'}, 'FontSize', 16, 'Location', 'northeast', 'Box', 'off')

text(16,0,  'Dark Box ''on''', 'VerticalAlignment', 'bottom', 'FontSize', 18)
ylim([0 100])

SetFigure([], gcf)


figure(5)
cla
hold on
rectangle('Position', [28, 0, 9, 100], 'FaceColor', [c_ord(7,:) .2], 'EdgeColor',[c_ord(7,:) .2])

h = errorbar(day(29:end), ArchT_vte_r_u(29:end), ArchT_vte_r_sem(29:end), 'color', green, 'LineWidth', 2);
h2 = plot(day(29:end), ArchT_vte_r_u(29:end), 's',  'MarkerSize', 14, 'MarkerFaceColor',  green, 'MarkerEdgeColor',  green);

hc = errorbar(day(29:end), ctrl_vte_r_u(29:end), ctrl_vte_r_sem(29:end), 'color', [.7 .7 .7], 'LineWidth', 2);
hc2 = plot(day(29:end), ctrl_vte_r_u(29:end), 'o',  'MarkerSize', 14, 'MarkerFaceColor',  [.7 .7 .7], 'MarkerEdgeColor',  [.7 .7 .7]);

ylabel('% of trials with VTEs')
xlabel('session')
xlim([28.5 35.5])
legend([h2, hc2], {'ArchT(n=4)', 'Control (n=3)'}, 'FontSize', 16, 'Location', 'northeast', 'Box', 'off')

text(29,0,  'Start Arm ''on''', 'VerticalAlignment', 'bottom',  'FontSize', 18)
ylim([0 100])

SetFigure([], gcf)


%%  Correct trials
green = [60 200 50]/255;
c_ord = linspecer(12); 
figure(12)
% title('VTE per session')
cla
hold on
rectangle('Position', [16, 0, 11, 100], 'FaceColor', [.6 .6 .6 .2], 'EdgeColor',[.6 .6 .6 .2])
rectangle('Position', [29, 0, 7, 100], 'FaceColor', [c_ord(7,:) .2], 'EdgeColor',[c_ord(7,:) .2])

h = errorbar(day, ArchT_c_u, ArchT_c_sem, 'color', green, 'LineWidth', 2);
h2 = plot(day, ArchT_c_u, 's',  'MarkerSize', 14, 'MarkerFaceColor',  green, 'MarkerEdgeColor',  green);

hc = errorbar(day, ctrl_c_u, ctrl_c_sem, 'color', [.7 .7 .7], 'LineWidth', 2);
hc2 = plot(day, ctrl_c_u, 'o',  'MarkerSize', 14, 'MarkerFaceColor',  [.7 .7 .7], 'MarkerEdgeColor',  [.7 .7 .7]);

ylabel('% correct')
xlabel('session')
xlim([0 35])
xline(15,'--k')
xline(28,'--k')
legend([h2, hc2], {'ArchT(n=4)', 'Control (n=3)'}, 'FontSize', 16, 'Location', 'southeast', 'Box', 'off')

text(0,50,  'Transparent Box ''on''', 'VerticalAlignment', 'bottom', 'FontSize', 18)
text(16,50,  'Dark Box ''on''', 'VerticalAlignment', 'bottom', 'FontSize', 18)
text(29,50,  'Start Arm ''on''', 'VerticalAlignment', 'bottom',  'FontSize', 18)
ylim([50 100])

SetFigure([], gcf)

figure(13)
cla
hold on

h = errorbar(day(3:14), ArchT_c_u(3:14), ArchT_c_sem(3:14), 'color', green, 'LineWidth', 2);
h2 = plot(day(3:14), ArchT_c_u(3:14), 's',  'MarkerSize', 14, 'MarkerFaceColor',  green, 'MarkerEdgeColor',  green);

hc = errorbar(day(3:14), ctrl_c_u(3:14), ctrl_c_sem(3:14), 'color', [.7 .7 .7], 'LineWidth', 2);
hc2 = plot(day(3:14), ctrl_c_u(3:14), 'o',  'MarkerSize', 14, 'MarkerFaceColor',  [.7 .7 .7], 'MarkerEdgeColor',  [.7 .7 .7]);

ylabel('% correct')
xlabel('session')
xlim([0 14.5])
legend([h2, hc2], {'ArchT(n=4)', 'Control (n=3)'}, 'FontSize', 16, 'Location', 'southeast', 'Box', 'off')
ylim([50 100])
text(0,50,  'Transparent Box ''on''', 'VerticalAlignment', 'bottom', 'FontSize', 18)

SetFigure([], gcf)

figure(14)
cla
hold on
rectangle('Position', [14, 0, 14, 100], 'FaceColor', [.6 .6 .6 .2], 'EdgeColor',[.6 .6 .6 .2])

h = errorbar(day(15:27), ArchT_c_u(15:27), ArchT_c_sem(15:27), 'color', green, 'LineWidth', 2);
h2 = plot(day(15:27), ArchT_c_u(15:27), 's',  'MarkerSize', 14, 'MarkerFaceColor',  green, 'MarkerEdgeColor',  green);

hc = errorbar(day(15:27), ctrl_c_u(15:27), ctrl_c_sem(15:27), 'color', [.7 .7 .7], 'LineWidth', 2);
hc2 = plot(day(15:27), ctrl_c_u(15:27), 'o',  'MarkerSize', 14, 'MarkerFaceColor',  [.7 .7 .7], 'MarkerEdgeColor',  [.7 .7 .7]);

ylabel('% correct')
xlabel('session')
xlim([15.5 27.5])
legend([h2, hc2], {'ArchT(n=4)', 'Control (n=3)'}, 'FontSize', 16, 'Location', 'southeast', 'Box', 'off')

text(16,50,  'Dark Box ''on''', 'VerticalAlignment', 'bottom', 'FontSize', 18)
ylim([50 100])

SetFigure([], gcf)


figure(15)
cla
hold on
rectangle('Position', [28, 0, 9, 100], 'FaceColor', [c_ord(7,:) .2], 'EdgeColor',[c_ord(7,:) .2])

h = errorbar(day(29:end), ArchT_c_u(29:end), ArchT_c_sem(29:end), 'color', green, 'LineWidth', 2);
h2 = plot(day(29:end), ArchT_c_u(29:end), 's',  'MarkerSize', 14, 'MarkerFaceColor',  green, 'MarkerEdgeColor',  green);

hc = errorbar(day(29:end), ctrl_c_u(29:end), ctrl_c_sem(29:end), 'color', [.7 .7 .7], 'LineWidth', 2);
hc2 = plot(day(29:end), ctrl_c_u(29:end), 'o',  'MarkerSize', 14, 'MarkerFaceColor',  [.7 .7 .7], 'MarkerEdgeColor',  [.7 .7 .7]);

ylabel('% correct')
xlabel('session')
xlim([28.5 35.5])
legend([h2, hc2], {'ArchT(n=4)', 'Control (n=3)'}, 'FontSize', 16, 'Location', 'southeast', 'Box', 'off')

text(29,50,  'Start Arm ''on''', 'VerticalAlignment', 'bottom',  'FontSize', 18)
ylim([50 100])

SetFigure([], gcf)
%% all subjects

figure(20)
% title('VTE per session')
cla
hold on
rectangle('Position', [16, 0, 11, 100], 'FaceColor', [.6 .6 .6 .2], 'EdgeColor',[.6 .6 .6 .2])
rectangle('Position', [29, 0, 7, 100], 'FaceColor', [c_ord(7,:) .2], 'EdgeColor',[c_ord(7,:) .2])

h = plot(day, ArchT_vte_r, 'color', green, 'LineWidth', 2);
h2 = plot(day, ArchT_vte_r, 's',  'MarkerSize', 14, 'MarkerFaceColor',  green, 'MarkerEdgeColor',  green);

hc =plot(day, ctrl_vte_r, 'color', [.7 .7 .7], 'LineWidth', 2);
hc2 = plot(day, ctrl_vte_r, 'o',  'MarkerSize', 14, 'MarkerFaceColor',  [.7 .7 .7], 'MarkerEdgeColor',  [.7 .7 .7]);

ylabel('% of trials with VTEs')
xlabel('session')
xlim([0 35])
xline(15,'--k')
xline(28,'--k')
legend([h2, hc2], {'ArchT(n=4)', 'Control (n=3)'}, 'FontSize', 16, 'Location', 'northeast', 'Box', 'off')

text(0,0,  'Transparent Box ''on''', 'VerticalAlignment', 'bottom', 'FontSize', 18)
text(16,0,  'Dark Box ''on''', 'VerticalAlignment', 'bottom', 'FontSize', 18)
text(29,0,  'Start Arm ''on''', 'VerticalAlignment', 'bottom',  'FontSize', 18)
ylim([0 100])

SetFigure([], gcf)



%% convert data to table
% subs = [ones
% VTE_t = table(



