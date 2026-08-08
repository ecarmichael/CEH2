%% Pox_HF_LFP

% loops over sessions and get the the bandpower values for bands of
% interest. 


fig_size = [100, 100, 600, 400]; % for consistent figure sizes. 

c_ord = MS_linspecer(3); 

%%

f_list = dir('pox*.mat'); 
for ii = length(f_list):-1:1
    s_list{ii} = f_list(ii).name(1:end-4); 
end
pox_idx = nan(size(s_list)); 
pow_mat_ca = NaN(length(s_list),5); 
pow_mat_sub = NaN(length(s_list),5); 

psd_ca = []; psd_sub = []; 


% index the pox mice(mouse)

for ii = length(s_list):-1:1

    load([s_list{ii} '.mat'])

    f = fieldnames(data); 
    this_sess = data.(f{1}); 

    if sum(contains(s_list{ii}, '3265')) > 0
        pox_idx(ii) = true;
    elseif sum(contains(s_list{ii}, '3265')) == 0
        pox_idx(ii) = false;
    end
   
    this_csc = this_sess.csc; 

    csc_mov = restrict(this_csc, this_sess.mov_iv); 

    csc_no_mov = restrict(this_csc, InvertIV(this_sess.mov_iv, this_csc.tvec(1), this_csc.tvec(end))); 


    %% grab the power from movement or non-movement epochs. 
    % overall
    pow_mat_ca(ii,1) = bandpower(this_csc.data(1,:), this_csc.cfg.hdr{1}.SamplingFrequency, [1 40]);
    %theta
    pow_mat_ca(ii,2) = bandpower(csc_mov.data(1,:), csc_mov.cfg.hdr{1}.SamplingFrequency, [6 10]); 
    % low gamma
    pow_mat_ca(ii,3) = bandpower(csc_mov.data(1,:), csc_mov.cfg.hdr{1}.SamplingFrequency, [30 45]); 
    %mid gamma
    pow_mat_ca(ii,4) = bandpower(csc_mov.data(1,:), csc_mov.cfg.hdr{1}.SamplingFrequency, [60 90]); 
    % ripples
    pow_mat_ca(ii,5) = bandpower(csc_no_mov.data(1,:), csc_no_mov.cfg.hdr{1}.SamplingFrequency, [120 200]); 


      % Sub overall
    pow_mat_sub(ii,1) = bandpower(this_csc.data(2,:), this_csc.cfg.hdr{1}.SamplingFrequency, [1 40]);
    %theta
    pow_mat_sub(ii,2) = bandpower(csc_mov.data(2,:), csc_mov.cfg.hdr{1}.SamplingFrequency, [6 10]); 
    % low gamma
    pow_mat_sub(ii,3) = bandpower(csc_mov.data(2,:), csc_mov.cfg.hdr{1}.SamplingFrequency, [30 45]); 
    %mid gamma
    pow_mat_sub(ii,4) = bandpower(csc_mov.data(2,:), csc_mov.cfg.hdr{1}.SamplingFrequency, [60 90]); 
    % ripples
    pow_mat_sub(ii,5) = bandpower(csc_no_mov.data(2,:), csc_no_mov.cfg.hdr{1}.SamplingFrequency, [120 200]); 

    %% grab a psd for plotting. 

    hann_win = 2^12; % always make this in base 2 for speed

    [psd_ca{ii}.pxx, psd_ca{ii}.f] = pwelch(this_csc.data(1,:), hanning(hann_win), hann_win/2, hann_win*2 , this_csc.cfg.hdr{1}.SamplingFrequency);

    [psd_sub{ii}.pxx, psd_sub{ii}.f] = pwelch(this_csc.data(2,:), hanning(hann_win), hann_win/2, hann_win*2 , this_csc.cfg.hdr{1}.SamplingFrequency);

end

pox_idx = logical(pox_idx); 

%%  quick Psd comparisons

figure(101)
set(gcf,'Units','pixels','position',fig_size); % helpful to keep the figures the same size. 

clf

subplot(2,3,1)
hold on

for ii = length(psd_ca):-1:1
    pxx_ca(ii,:) = 10*log10(psd_ca{ii}.pxx)./max(10*log10(psd_ca{ii}.pxx)); 
end


% remove the notch indices for a cleaner plot. 
notch_idx = nearest_idx([58 62 118 122 178 182 238 242], psd_ca{ii}.f);
keep_idx = ones(size(psd_ca{ii}.f)); 
keep_idx([notch_idx(1):notch_idx(2), notch_idx(3):notch_idx(4), notch_idx(5):notch_idx(6), notch_idx(7):notch_idx(8)]) = 0;
keep_idx = logical(keep_idx); 

% err = (std(pxx(~pox_idx,:),0,1, 'omitmissing')./sqrt(size(pxx(~pox_idx,:),1)))
h = shadedErrorBar(psd_ca{ii}.f(keep_idx), mean(pxx_ca(~pox_idx,keep_idx)),MS_SEM_vec(pxx_ca(~pox_idx,keep_idx)));
h.mainLine.LineWidth = 2; 

% plot the Pox
h = shadedErrorBar(psd_ca{ii}.f(keep_idx), mean(pxx_ca(pox_idx,keep_idx)),MS_SEM_vec(pxx_ca(pox_idx,keep_idx))); 
h.mainLine.Color = c_ord(2,:); 
h.mainLine.LineWidth = 2; 
h.patch.FaceColor = h.mainLine.Color; 
h.patch.EdgeColor = h.mainLine.Color;
h.edge(1).Color = h.mainLine.Color; 
h.edge(2).Color = h.mainLine.Color; 

% Set axis labels and title for the plot
xlabel('Frequency (Hz)');
ylabel('Normalized Power');
title('Ca1 PSD Comparison');
legend({'Control', 'Pox'}, 'Location', 'Best', 'Box','off');

xlim([0 250])

set(gca, 'fontsize', 8)

% same for the Subiculum


subplot(2,3,4)
hold on

for ii = length(psd_sub):-1:1
    pxx_sub(ii,:) = 10*log10(psd_sub{ii}.pxx)./max(10*log10(psd_sub{ii}.pxx)); 
end


% remove the notch indices for a cleaner plot. 
notch_idx = nearest_idx([58 62 118 122 178 182 238 242], psd_sub{ii}.f);
keep_idx = ones(size(psd_sub{ii}.f)); 
keep_idx([notch_idx(1):notch_idx(2), notch_idx(3):notch_idx(4), notch_idx(5):notch_idx(6), notch_idx(7):notch_idx(8)]) = 0;
keep_idx = logical(keep_idx); 

% err = (std(pxx(~pox_idx,:),0,1, 'omitmissing')./sqrt(size(pxx(~pox_idx,:),1)))
h = shadedErrorBar(psd_sub{ii}.f(keep_idx), mean(pxx_sub(~pox_idx,keep_idx)),MS_SEM_vec(pxx_sub(~pox_idx,keep_idx)));
h.mainLine.LineWidth = 2; 

% plot the Pox
h = shadedErrorBar(psd_sub{ii}.f(keep_idx), mean(pxx_sub(pox_idx,keep_idx)),MS_SEM_vec(pxx_sub(pox_idx,keep_idx))); 
h.mainLine.Color = c_ord(2,:); 
h.mainLine.LineWidth = 2; 
h.patch.FaceColor = h.mainLine.Color; 
h.patch.EdgeColor = h.mainLine.Color;
h.edge(1).Color = h.mainLine.Color; 
h.edge(2).Color = h.mainLine.Color; 

% Set axis labels and title for the plot
xlabel('Frequency (Hz)');
ylabel('Normalized Power');
title('Subiculum PSD Comparison');
legend({'Control', 'Pox'}, 'Location', 'Best', 'Box','off');

xlim([0 250])

set(gca, 'fontsize', 8)

%% plot the power values
subplot(2,3,2); 
cla
hold on
x_mat = [1, 2; 3, 4; 5, 6; 7, 8]; 
% plot them
for ii = 2:5
    MS_bar_w_err(pow_mat_ca(~pox_idx,ii)./pow_mat_ca(~pox_idx,1), pow_mat_ca(pox_idx,ii)./pow_mat_ca(pox_idx,1), [.8 .8 .8; c_ord(2,:)], 1, 'ttest2', x_mat(ii-1,:))
end

% Finalize the power plot with appropriate labels and settings
xlabel('Power Bands');
ylabel('Normalized Power');
title('Power Values Comparison');
set(gca, 'XTick', 1.5:2:8.5,'TickLabelInterpreter','tex', 'XTickLabel', { 'Theta', 'Low\newlineGamma',  'Mid\newlineGamma', 'Ripples'});
set(gca, 'fontsize', 8);


subplot(2,3,3); 

% anovan 
pow_mat_norm = pow_mat_ca(:,2:end)./pow_mat_ca(:,1); 
type = categorical([ones(1,size(pow_mat_norm,1)), 2*ones(1,size(pow_mat_norm,1)), 3*ones(1,size(pow_mat_norm,1)), 4*ones(1,size(pow_mat_norm,1))])';
pox = logical([pox_idx, pox_idx, pox_idx, pox_idx]); 

[~,~,stats] = anovan([pow_mat_norm(:,1); pow_mat_norm(:,2); pow_mat_norm(:,3); pow_mat_norm(:,4)],...
    {type pox},'model',2,'varnames',{'band','pox'}); 

[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A")=gnames(tbl.("Group A"));
tbl.("Group B")=gnames(tbl.("Group B"));


subplot(2,3,5); 
cla
hold on
x_mat = [1, 2; 3, 4; 5, 6; 7, 8]; 
% plot them
for ii = 2:5
    MS_bar_w_err(pow_mat_sub(~pox_idx,ii)./pow_mat_sub(~pox_idx,1), pow_mat_sub(pox_idx,ii)./pow_mat_sub(pox_idx,1), [.8 .8 .8; c_ord(2,:)], 1, 'ttest2', x_mat(ii-1,:))
end

% Finalize the power plot with appropriate labels and settings
xlabel('Power Bands');
ylabel('Normalized Power');
title('Subiculum Power Comparison');
set(gca, 'XTick', 1.5:2:8.5,'TickLabelInterpreter','tex', 'XTickLabel', { 'Theta', 'Low\newlineGamma',  'Mid\newlineGamma', 'Ripples'});
set(gca, 'fontsize', 8);


subplot(2,3,6); 

% anovan 
pow_mat_norm = pow_mat_sub(:,2:end)./pow_mat_sub(:,1); 
type = categorical([ones(1,size(pow_mat_norm,1)), 2*ones(1,size(pow_mat_norm,1)), 3*ones(1,size(pow_mat_norm,1)), 4*ones(1,size(pow_mat_norm,1))])';
pox = logical([pox_idx, pox_idx, pox_idx, pox_idx]); 

[~,~,stats] = anovan([pow_mat_norm(:,1); pow_mat_norm(:,2); pow_mat_norm(:,3); pow_mat_norm(:,4)],...
    {type pox},'model',2,'varnames',{'band','pox'}); 

[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A")=gnames(tbl.("Group A"));
tbl.("Group B")=gnames(tbl.("Group B"));

%% same thing but for Ca1
clear

load('all_data.mat')

% keep only the TL sessions

s_list = fieldnames(all_data);
% only keep the LT or TFC days
k_idx = false(size(s_list));
for ii = 1:length(s_list)
    if contains(s_list{ii}, 'tl')
        k_idx(ii) = true;
    end
    if contains(s_list{ii}, '2217_tl4')
        k_idx(ii) = false;
    end
end

s_list(~k_idx) = []; 
pox_idx = nan(size(s_list)); 
pow_mat_ca = NaN(length(s_list),5); 
psd = []; 

% index the pox mice(mouse)

for ii = length(s_list):-1:1

    if sum(contains(s_list{ii}, '3265')) > 0
        pox_idx(ii) = true;
    elseif sum(contains(s_list{ii}, '3265')) == 0
        pox_idx(ii) = false;
    end
   
    this_csc = all_data.(s_list{ii}).csc; 

    % get the movement from the encoder

    move_ts = ts({all_data.(s_list{ii}).evts.t{contains(all_data.(s_list{ii}).evts.label, '8')}});

    mov_rate  = MS_spike2rate(move_ts, this_csc.tvec);

    cfg_mov = [];
    cfg_mov.threshold = .001;
    cfg_mov.dcn = '<';
    cfg_mov.operation = '<';
    cfg_mov.minlen = .5;

    mov_iv = TSDtoIV(cfg_mov, mov_rate);

    % pad movement

    cfg_resize.d = [+.5 -.5];

    mov_iv = ResizeIV(cfg_resize, mov_iv);

    csc_mov = restrict(this_csc, mov_iv); 

    csc_no_mov = restrict(this_csc, InvertIV(mov_iv, this_csc.tvec(1), this_csc.tvec(end))); 


    %% grab the power from movement or non-movement epochs. 
    % overall
    pow_mat_ca(ii,1) = bandpower(this_csc.data, this_csc.cfg.hdr{1}.SamplingFrequency, [1 40]);
    %theta
    pow_mat_ca(ii,2) = bandpower(csc_mov.data, csc_mov.cfg.hdr{1}.SamplingFrequency, [6 10]); 
    % low gamma
    pow_mat_ca(ii,3) = bandpower(csc_mov.data, csc_mov.cfg.hdr{1}.SamplingFrequency, [30 45]); 
    %mid gamma
    pow_mat_ca(ii,4) = bandpower(csc_mov.data, csc_mov.cfg.hdr{1}.SamplingFrequency, [60 90]); 
    % ripples
    pow_mat_ca(ii,5) = bandpower(csc_no_mov.data, csc_no_mov.cfg.hdr{1}.SamplingFrequency, [120 200]); 

    %% grab a psd for plotting. 

    hann_win = 2^12; % always make this in base 2 for speed

    [psd{ii}.pxx, psd{ii}.f] = pwelch(this_csc.data, hanning(hann_win), hann_win/2, hann_win*2 , this_csc.cfg.hdr{1}.SamplingFrequency);


end

pox_idx = logical(pox_idx); 

%%  quick Psd comparisons

figure(100)
set(gcf,'Units','pixels','position',fig_size); % helpful to keep the figures the same size. 

clf

subplot(2,2,1)
hold on
ctrl_pxx = []; pox_pxx = []; 

for ii = length(psd):-1:1
    pxx(ii,:) = 10*log10(psd{ii}.pxx)./max(10*log10(psd{ii}.pxx)); 
end



% plot(psd{ii}.f, mean(pxx(~pox_idx,:)), 'k', 'LineWidth',2)
% plot(psd{ii}.f, mean(pxx(pox_idx,:)), 'r', 'LineWidth',2)

% remove the notch indices for a cleaner plot. 
notch_idx = nearest_idx([58 62 118 122 178 182 238 242], psd{ii}.f);
keep_idx = ones(size(psd{ii}.f)); 
keep_idx([notch_idx(1):notch_idx(2), notch_idx(3):notch_idx(4), notch_idx(5):notch_idx(6), notch_idx(7):notch_idx(8)]) = 0;
keep_idx = logical(keep_idx); 

% err = (std(pxx(~pox_idx,:),0,1, 'omitmissing')./sqrt(size(pxx(~pox_idx,:),1)))
h = shadedErrorBar(psd{ii}.f(keep_idx), mean(pxx(~pox_idx,keep_idx)),MS_SEM_vec(pxx(~pox_idx,keep_idx)));
h.mainLine.LineWidth = 2; 

% plot the Pox
h = shadedErrorBar(psd{ii}.f(keep_idx), mean(pxx(pox_idx,keep_idx)),MS_SEM_vec(pxx(pox_idx,keep_idx))); 
h.mainLine.Color = c_ord(2,:); 
h.mainLine.LineWidth = 2; 
h.patch.FaceColor = h.mainLine.Color; 
h.patch.EdgeColor = h.mainLine.Color;
h.edge(1).Color = h.mainLine.Color; 
h.edge(2).Color = h.mainLine.Color; 



% Set axis labels and title for the plot
xlabel('Frequency (Hz)');
ylabel('Normalized Power');
title('Power Spectral Density Comparison');
legend({'Control', 'Pox'}, 'Location', 'Best', 'Box','off');

xlim([0 250])

set(gca, 'fontsize', 8)

title('CA1')

% plot the power values
subplot(2,2,2); 
cla
hold on
x_mat = [1, 2; 3, 4; 5, 6; 7, 8]; 
% plot them
for ii = 2:5
    MS_bar_w_err(pow_mat_ca(~pox_idx,ii)./pow_mat_ca(~pox_idx,1), pow_mat_ca(pox_idx,ii)./pow_mat_ca(pox_idx,1), [.8 .8 .8; c_ord(2,:)], 1, 'ttest2', x_mat(ii-1,:))
end

% Finalize the power plot with appropriate labels and settings
xlabel('Power Bands');
ylabel('Normalized Power');
title('Power Values Comparison');
set(gca, 'XTick', 1.5:2:8.5,'TickLabelInterpreter','tex', 'XTickLabel', { 'Theta', 'Low\newlineGamma',  'Mid\newlineGamma', 'Ripples'});
set(gca, 'fontsize', 8);


subplot(2,2,4); 

% anovan 
pow_mat_norm = pow_mat_ca(:,2:end)./pow_mat_ca(:,1); 
type = categorical([ones(1,size(pow_mat_norm,1)), 2*ones(1,size(pow_mat_norm,1)), 3*ones(1,size(pow_mat_norm,1)), 4*ones(1,size(pow_mat_norm,1))])';
pox = logical([pox_idx; pox_idx; pox_idx; pox_idx]); 

[p,~,stats] = anovan([pow_mat_norm(:,1); pow_mat_norm(:,2); pow_mat_norm(:,3); pow_mat_norm(:,4)],...
    {type pox},'model',2,'varnames',{'band','pox'}); 

[results,~,~,gnames] = multcompare(stats,"Dimension",[1 2]);

tbl = array2table(results,"VariableNames", ...
    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);
tbl.("Group A")=gnames(tbl.("Group A"));
tbl.("Group B")=gnames(tbl.("Group B"));
