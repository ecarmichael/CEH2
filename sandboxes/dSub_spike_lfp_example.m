%% example of spikes with LFP

cd('/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/dSubiculum/Incoming/M23/M23_2021-08-08_D13')

Meta = MS_Load_meta; 
%%
cfg.fc = {'TT4_01_5.t', 'TT4_03_2.t', 'TT4_02_4.t'};
cfg.getTTnumbers = 0; 

S = LoadSpikes(cfg);

cfg = [];
cfg.fc = {Meta.goodCSC};

csc = LoadCSC(cfg);


evt = LoadEvents([]);


csc_r  =  restrict(csc, evt.t{2}(2), evt.t{3}(2));

S_r  =  restrict(S, evt.t{2}(2), evt.t{3}(2));

%% 
close all
figure(101)
t1 = 60*2000; %97598
t2 = 110*2000; %99598

% t1 = 97598+800;
% t2 = 99598;

plot(S_r)
hold on

plot(csc_r.tvec, csc_r.data*5000, 'k')
xlim([csc_r.tvec(t1) csc_r.tvec(t2)]); 


%% restrict to a specific window

S_r2 = restrict(S_r, csc_r.tvec(t1), csc_r.tvec(t2));
% 
csc_r2 = restrict(csc_r, csc_r.tvec(t1), csc_r.tvec(t2));


cfg_theta = [];
cfg_theta.f = [6 12];
theta = FilterLFP(cfg_theta, csc_r2);

% S_r2 = S_r;
% csc_r2 = csc_r;
%% plot again but with spikes on the LFP

c_ord = linspecer(length(S.t)+2); 
close all
figure(120)

c1_idx = nearest_idx(S_r2.t{1}, csc_r2.tvec);
c2_idx = nearest_idx(S_r2.t{2}, csc_r2.tvec);
c3_idx = nearest_idx(S_r2.t{3}, csc_r2.tvec);

plot(csc_r2.tvec, csc_r2.data, '-k', 'linewidth',1.5)
hold on
plot(theta.tvec, theta.data - 0.0005, 'color', c_ord(5,:), 'linewidth',3) 

hold on
for ii = 1:length(c1_idx)
    plot([csc_r2.tvec(c1_idx(ii))' csc_r2.tvec(c1_idx(ii))'], [csc_r2.data(c1_idx(ii)), csc_r2.data(c1_idx(ii))+.00025], 'color', c_ord(1,:), 'linewidth',1.5)
    plot(csc_r2.tvec(c1_idx(ii)), csc_r2.data(c1_idx(ii))+.00026,'o', 'color', c_ord(1,:), 'linewidth',3)
    
    plot([theta.tvec(c1_idx(ii))' theta.tvec(c1_idx(ii))'], [theta.data(c1_idx(ii)) - 0.0005, theta.data(c1_idx(ii))+.00025 - 0.0005], 'color', c_ord(1,:), 'linewidth',3)
end

for ii = 1:length(c2_idx)
    plot([csc_r2.tvec(c2_idx(ii))' csc_r2.tvec(c2_idx(ii))'], [csc_r2.data(c2_idx(ii)), csc_r2.data(c2_idx(ii))+.00025], 'color', c_ord(2,:), 'linewidth',1.5)
    plot(csc_r2.tvec(c2_idx(ii)), csc_r2.data(c2_idx(ii))+.00026,'o', 'color', c_ord(2,:), 'linewidth',3)    
    
    plot([theta.tvec(c2_idx(ii))' theta.tvec(c2_idx(ii))'], [theta.data(c2_idx(ii)) - 0.0005, theta.data(c2_idx(ii))+.00025 - 0.0005], 'color', c_ord(2,:), 'linewidth',3)
end

for ii = 1:length(c3_idx)
    plot([csc_r2.tvec(c3_idx(ii))' csc_r2.tvec(c3_idx(ii))'], [csc_r2.data(c3_idx(ii)), csc_r2.data(c3_idx(ii))+.00025], 'color', c_ord(3,:), 'linewidth',1.5)
    plot(csc_r2.tvec(c3_idx(ii)), csc_r2.data(c3_idx(ii))+.00026,'o', 'color', c_ord(3,:), 'linewidth',3)    
end

plot(csc_r2.tvec, csc_r2.data, '-k', 'linewidth',1.5)

% xlim([8.8274e+03 (8.8274e+03)+.5])
xlim([8.8018e+03 (8.8018e+03)+.25])