%% sandbox_NAMI_SWR_analyses 



% Initialize some things: 


c_ord = MS_linspecer(5); %generate a few nice colours. 

fig_size = [100, 100, 500, 400]; % for consistent figure sizes. 
ft_size = 12;  %font size

%% loop over all the sessions of interest. Ultimately you will be able to loop over the sessions and hold the metrics of interest in order to get some session/subject level comparisons. 

% but this is for later. 

%% start with a nice session and generate some simple plots.  
fname = 'pox3568_TFCD5'; % 3568 TFCD2/4/5 are all nice. 

this_sess = load([fname '.mat']); 

this_sess = this_sess.data.(fname); % get the nice data from this session 

%% Get the relative timing of each Sub SWR compared to the nearest Ca1 SWR

% get the center of the events
sub_cent = IVcenters(this_sess.swrs_sub);
ca1_cent = IVcenters(this_sess.swrs_ca1);

% loop over the events and find pairs and their offest

pair_swr = []; % make an empty place to put things. 
offset_max = .5; % max time between CA1 peak and Sub peak + or minus. 

for ii = length(sub_cent):-1:1
    this_d = sub_cent(ii) - ca1_cent;
    co_evt = find(abs(this_d) < .50);

    if isempty(co_evt) || length(co_evt) > 1
        pair_swr(ii) = NaN;
    else
        pair_swr(ii) = this_d(co_evt);
    end
end


% make a new set of events for Ca1 leading or Sub leading 
this_sess.swrs_sub_atyp = SelectIV([], this_sess.swrs_sub, pair_swr < 0); 
this_sess.swrs_sub_std = SelectIV([], this_sess.swrs_sub, pair_swr > 0); 

% find the CA1 SWRs that overlap with each Sub SWR type
this_sess.swrs_ca1_atyp  = IntersectIV([], this_sess.swrs_ca1, this_sess.swrs_sub_atyp); 
this_sess.swrs_ca1_std  = IntersectIV([], this_sess.swrs_ca1, this_sess.swrs_sub_std); 

%% make an event triggered spectrogram of the CA1 and Sub SWRs relative to the CA1 center using fieldtrip. 

% skip this if you don't want to deal with fieldtrip.

% get things ready for field trip.  You will need to clone the fieldtrip
% repo and initialize it by moving to the fieldtrip folder in matlab and
% running ft_defaults. 

% pull Ca1 channel out first. 

csc_ft = this_sess.csc;
csc_ft.data = this_sess.csc.data(1,:);
csc_ft.label = [];
csc_ft.label{1} = this_sess.csc.label{1};

[csc_ft_out, TFR] = Triggered_Spec_FT(csc_ft, ca1_cent, [], 80:.5:200, [], [-.25 .25], 1);
set(gcf,'Units','pixels','position',fig_size); % helpful to keep the figures the same size. 

xlim([-.05 .05]); ylabel('frequency (hz)'); xlabel('time from CA1 SWR center (ms)')
set(gca,'xtick', -.05:0.025:.05);
set(gca, 'XTickLabel', get(gca, 'XTick')*1000)
title('Ca1 SWRS')
axis("square")

% move this to a subplot. 
% ax = gca;
% figure(3000)
% subplot(3,2,1, ax)

% Next Sub aligned to CA1 

csc_ft = this_sess.csc;
csc_ft.data = this_sess.csc.data(2,:);
csc_ft.label = [];
csc_ft.label{1} = this_sess.csc.label{2};

[csc_ft_out, TFR] = Triggered_Spec_FT(csc_ft, ca1_cent, [], 80:.5:200, [], [-.25 .25], 1);
set(gcf,'Units','pixels','position',fig_size); % helpful to keep the figures the same size. 

xlim([-.05 .05]); ylabel('frequency (hz)'); xlabel('time from CA1 SWR center (ms)')
set(gca,'xtick', -.05:0.025:.05);
set(gca, 'XTickLabel', get(gca, 'XTick')*1000)
title('Sub SWRS')
axis("square")


% Next Sub LEF aligned to CA1 standard SWR
ca1_std_cent = IVcenters(this_sess.swrs_ca1_std); 

[csc_ft_out, TFR] = Triggered_Spec_FT(csc_ft, ca1_std_cent, [], 80:.5:200, [], [-.25 .25], 1);
set(gcf,'Units','pixels','position',fig_size); % helpful to keep the figures the same size. 

xlim([-.05 .05]); ylabel('frequency (hz)'); xlabel('time from CA1 SWR center (ms)')
set(gca,'xtick', -.05:0.025:.05);
set(gca, 'XTickLabel', get(gca, 'XTick')*1000)
title('Sub SWRS - Standard')
axis("square");


% Next Sub LEF aligned to CA1 atypical SWR
ca1_atyp_cent = IVcenters(this_sess.swrs_ca1_atyp); 

[csc_ft_out, TFR] = Triggered_Spec_FT(csc_ft, ca1_atyp_cent, [], 80:.5:200, [], [-.25 .25], 1);
set(gcf,'Units','pixels','position',fig_size); % helpful to keep the figures the same size. 

xlim([-.05 .05]); ylabel('frequency (hz)'); xlabel('time from CA1 SWR center (ms)')
set(gca,'xtick', -.05:0.025:.05);
set(gca, 'XTickLabel', get(gca, 'XTick')*1000)
title('Sub SWRS - Atypical')
axis("square")



%% make a figure showing the distribution of CA1-Sub SWR timing within a session. 
% subplot(3,2,5) % subplot 
figure(3002); clf
set(gcf,'Units','pixels','position',fig_size); % helpful to keep the figures the same size. 

hold on
all_pairs = pair_swr(~isnan(pair_swr))*1000; %take the non-nan pair times and converting them to miliseconds

pos_pairs = all_pairs(all_pairs> 0); 
neg_pairs = all_pairs(all_pairs< 0); 

% make two histograms (one for positive and one for negative so that they
% can have different colours

histogram(pos_pairs, -25:1:25, 'FaceColor',c_ord(1,:), 'EdgeColor',c_ord(1,:)) % make a histogram  within +/- 25ms
histogram(neg_pairs, -25:1:25, 'FaceColor',c_ord(2,:), 'EdgeColor',c_ord(2,:)) % make a histogram  within +/- 25ms

xline(0, '--k', 'LineWidth',1)
% axis('square') % make things square. 

% Calculate and display the mean and standard deviation of the offsets
meanOffset = mean(pair_swr, 'omitmissing');
stdOffset = std(pair_swr, 'omitmissing');
title(sprintf('Mean: %.2f ms, Std: %.2f ms', meanOffset * 1000, stdOffset * 1000), 'FontSize',ft_size, 'FontWeight','normal');

% set the font to be a consistent size. 

set(gca, 'FontSize',ft_size)
axis("square")


%% Quantify the occurence of both types of SWRS within the task phases. 
% same thing you had done previsously using the restrict function, but this
% time instead of printing the result you will want to hold the values in a
% variable for later comparison.  YOu will need to add in the standard ad
% atypical swrs too. 

% I would do this an array with some structure like this (rows are
% conditions (pre, post, reward,...) and columns are the 4 types of SWR.  

swr_counts = []; % empty it to start. 
swr_counts_labels = {'Ca1', 'Sub', 'Std', 'Atyp'}; % labels for later. 

% restrict data to the pre task baseline
% ca1
swr_ca1_pre = restrict(this_sess.swrs_ca1, this_sess.csc.tvec(1), this_sess.evts.t{2}(1)); 
%sub 
swr_sub_pre = restrict(this_sess.swrs_sub, this_sess.csc.tvec(1), this_sess.evts.t{2}(1)); 
swr_sub_std_pre = restrict(this_sess.swrs_sub_std, this_sess.csc.tvec(1), this_sess.evts.t{2}(1)); 
swr_sub_atyp_pre = restrict(this_sess.swrs_sub_atyp, this_sess.csc.tvec(1), this_sess.evts.t{2}(1)); 


fprintf('Ca1 SWRs in Pre: <strong>%d</strong> \n', length(swr_ca1_pre.tstart))
fprintf('Sub SWRs in Pre: <strong>%d</strong> \n', length(swr_sub_pre.tstart))
fprintf('std Sub SWRs in Pre: <strong>%d</strong> \n', length(swr_sub_std_pre.tstart))
fprintf('atyp Sub SWRs in Pre: <strong>%d</strong> \n', length(swr_sub_atyp_pre.tstart))

% fill in the table. Coule be done in one line using a concatenation, but
% this keeps it easy. 
swr_counts(1,1) = length(swr_ca1_pre.tstart); 
swr_counts(1,2) = length(swr_sub_pre.tstart); 
swr_counts(1,3) = length(swr_sub_std_pre.tstart); 
swr_counts(1,4) = length(swr_sub_atyp_pre.tstart); 
% add a label for the rows
swr_counts_rows{1} = 'Pre'; 



% NAMI TO DO:  same as above for the other measures. 
% restrict data to the post period. 

% ca1
swr_ca1_post = restrict(this_sess.swrs_ca1, this_sess.evts.t{2}(2), this_sess.csc.tvec(end)); 
%sub 
swr_sub_post = restrict(this_sess.swrs_sub, this_sess.evts.t{2}(2), this_sess.csc.tvec(end)); 

fprintf('Ca1 SWRs in Post: <strong>%d</strong> \n', length(swr_ca1_post.tstart))

fprintf('Sub SWRs in Post: <strong>%d</strong> \n', length(swr_sub_post.tstart))


% get the baseline periods (30s before any tone) inside the TFC
dur = 30; 
tone_t_on = sort([this_sess.evts.t{4}(1:2:end); this_sess.evts.t{5}(1:2:end)]);

swr_ca1_base = restrict(this_sess.swrs_ca1, tone_t_on-30, tone_t_on); 
swr_sub_base = restrict(this_sess.swrs_sub, tone_t_on-30, tone_t_on); 

fprintf('Ca1 SWRs in baseline: <strong>%d</strong>   |    Sub SWRs in baseline: <strong>%d</strong> \n',...
    length(swr_ca1_base.tstart),  length(swr_sub_base.tstart))



% NAMI To Do:  same as above for the specific behaviour periods. 
% get the tone 1 periods 




% get the tone 2 periods



% tone 1 trace period (20 sec after end of tone)



% tone 2 trace period 