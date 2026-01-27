function [data_out] = HF_metrics(data_in, plot_flag, TTL)
%% HF_Metrics:  computes simple spiking measures

if nargin < 2
    plot_flag = 0;
    TTL = {'6', '7'};
elseif nargin < 3
    TTL = {'6', '7'};
end

% ToDo:
% - extract waveforms from phy2 wrapper.
% - classify cells using waveforms.


%% loop over cells
fr = []; isi =[]; b_idx = []; opto_resp = [];
for iS = length(data_in.S.t):-1:1
    S_metrics{iS} = [];
    S_metrics{iS}.fr = length(data_in.S.t{iS})./(data_in.csc.tvec(end)-data_in.csc.tvec(1));
    S_metrics{iS}.ISI = mean(diff(data_in.S.t{iS}));
    S_metrics{iS}.burst_idx = sum(diff(data_in.S.t{iS}) < .010) / sum(diff(data_in.S.t{iS}) > .010);

        
    % autocorrelation
    xbin_centers = -.025-0.001:0.001:0.025+0.001; % first and last bins are to be deleted later
    ac = zeros(size(xbin_centers));

    for iSpk = 1:length(data_in.S.t{iS})
        relative_spk_t = data_in.S.t{iS} - data_in.S.t{iS}(iSpk);
        ac = ac + hist(relative_spk_t,xbin_centers); % note that hist() puts all spikes outside the bin centers in the first and last bins! delete later.
    end

    xbin = xbin_centers(2:end-1); % remove unwanted bins
    zero_idx = xbin == 0;
    ac = ac(2:end-1);
    ac(zero_idx) = 0;

    S_metrics{iS}.auto_corr = ac;
    S_metrics{iS}.auto_corr_xbin = xbin*1000;

    % for internal plotting checks
    fr(iS) = S_metrics{iS}.fr;
    isi(iS) = S_metrics{iS}.ISI;
    b_idx(iS) = S_metrics{iS}.burst_idx;

    %% opto response
    c_red = [0.9153    0.2816    0.2878;
        0 0 1];
    
    TTL_name = {'red', 'blue'};

    for iO = 1:2
        if sum(contains(data_in.evts.label, TTL{iO}))
            this_TTL_idx = find(contains(data_in.evts.label, TTL{iO}));

            evt_t = data_in.evts.t{this_TTL_idx} ;

            % isolate events of a certain length
            e_d = evt_t(2,:) - evt_t(1,:);
            ITIs = unique(round(e_d, 3));
            for iTi = length(ITIs):-1:1

                if plot_flag
                    figure((iS*1000) + iO*10 + iTi)
                end
                evt_t = data_in.evts.t{this_TTL_idx} ;
                k_idx = round(e_d,3) == ITIs(iTi) ;
                evt_t(:,~k_idx) = [];


                if S_metrics{iS}.fr < 0.5
                    S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.dur = mode(e_d(k_idx));
                    S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.ITI = mode(diff(evt_t));
                    S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.pre = NaN;
                    S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.post = NaN;

                    S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.shuff_pre = NaN;
                    S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.shuff_post = NaN;

                    % if ~isempty(pre) || ~isempty(post)
                    S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.pval = NaN;

                    if iTi == 1
                        opto_resp.(TTL_name{iO})(iS) = NaN;
                    end
                end

                % isolate the cell
                S = data_in.S;
                S.t = [];
                S.t{1} = data_in.S.t{iS};
                S.label = [];
                S.label{1} =  data_in.S.label{iS};

                evt_t = sort(evt_t(1,:));

                cfg_peth = [];
                cfg_peth.window = [-.25 .25];
                cfg_peth. plot_type = 'raw';
                cfg_peth.dt = 0.001;
                cfg_peth.gauss_window = .025;
                cfg_peth.gauss_sd = .0025;
                cfg_peth.shuff = 500;
                cfg_peth.t_on = mode(e_d(k_idx));
                cfg_peth.rec_color = c_red(iO,:);
                if plot_flag
                    cfg_peth.plot = 'on';
                else
                    cfg_peth.plot = 'off';
                end
                [peth_S, ~,~, ~, ~, ~, ~, ~, ~,peth_T] = SpikePETH_Shuff(cfg_peth, S, evt_t);
                disp(S.label{1})
                % pre_fr = mean(zval(peth_IT<0));
                % laser_fr = mean(zval(find(peth_IT==0):nearest_idx(mode(e_d(k_idx)), peth_IT)));

                % t-test for pre and post rate for each stim?
                % u_val = unique(peth_T); pre = NaN(1,length(u_val)); post = pre;
                % for iV = length(u_val):-1:1
                %     this_idx = peth_T == u_val(iV);
                %     if isempty(this_idx)
                %         continue
                %     end
                %     pre(iV) = sum((peth_S(this_idx) > -mode(e_d(k_idx))) & (peth_S(this_idx) < 0));
                %     post(iV) = sum((peth_S(this_idx) > 0) & peth_S(this_idx) < mode(e_d(k_idx)));
                % end

                %rate based FR Stim changes
                this_fr = MS_spike2rate(S, data_in.csc.tvec, 0.001, 0.005); % needs small gau sd for precision.

                for ii = length(evt_t):-1:1
                    l_iv = iv(evt_t(ii), evt_t(ii)+ cfg_peth.t_on );
                    fr_r = restrict(this_fr, l_iv);
                    l_fr(ii) = mean(fr_r.data);

                    p_iv = iv(evt_t(ii)- cfg_peth.t_on, evt_t(ii) );
                    fr_r = restrict(this_fr, p_iv);
                    p_fr(ii) = mean(fr_r.data);

                end

                if plot_flag
                    figure(iS*1000+iO*10+iTi+100000)
                    MS_bar_w_err(p_fr, l_fr, [.6 .6 .6; cfg_peth.rec_color], 1, 'ttest');
                    set(gca, 'XTickLabel', {'baseline', 'light'})
                end
                % shuffles
                real_iv = iv(evt_t- cfg_peth.t_on , evt_t+ cfg_peth.t_on);
                shuff_t = datasample(this_fr.tvec, 10000, 1, "Replace",true);
                shuff_iv = iv(shuff_t, shuff_t+ cfg_peth.t_on);

                cfg_diff.verbose = 0;
                shuff_iv = DifferenceIV(cfg_diff, shuff_iv, real_iv);

                shuff_id = datasample(1:length(shuff_iv.tstart), 1000);
                shuff_l_fr = NaN(1000,1); shuff_p_fr = shuff_l_fr;

                for ishuff = 1:100
                    this_s_fr = restrict(this_fr, shuff_iv.tstart(shuff_id(ishuff)), shuff_iv.tend(shuff_id(ishuff)));
                    shuff_l_fr(ishuff) = mean(this_s_fr.data);

                    this_s_fr = restrict(this_fr, shuff_iv.tstart(shuff_id(ishuff)) - cfg_peth.t_on , shuff_iv.tstart(shuff_id(ishuff)));
                    shuff_p_fr(ishuff) = mean(this_s_fr.data);
                end


                S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.dur = mode(e_d(k_idx));
                S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.ITI = mode(diff(evt_t));
                S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.pre = p_fr;
                S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.post = l_fr;

                S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.shuff_pre = shuff_l_fr;
                S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.shuff_post = shuff_p_fr;

                % if ~isempty(pre) || ~isempty(post)
                [~,S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.pval] = ttest(l_fr, p_fr);

                if iTi == 1
                    opto_resp.(TTL_name{iO})(iS) = ttest(l_fr, p_fr);
                end
                % else
                % S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.pval = NaN;
                % opto_resp.(TTL_name{iO})(iS) = NaN;
            end
        end % end ITIS
    end
end% opto colours


%% collect the output stats

data_out = data_in; 
data_out.S_metrics  = S_metrics; 

%% quick check
r_resp = NaN(length(S_metrics), 1); 
for ii = 1:length(S_metrics)

    fprintf('Red Response Cell #%.d, %.4f\n',ii,  S_metrics{ii}.opto_red{1}.pval)
    r_resp(ii) = (S_metrics{ii}.opto_red{1}.pval < 0.05);
end

for iO = 1:length(fieldnames(opto_resp))
    fprintf('%s optoresponsive cells %.2f%% (%.d/%.d)\n',TTL_name{iO}, (sum(opto_resp.(TTL_name{iO}), 'omitmissing')./length(opto_resp.(TTL_name{iO})))*100, sum(opto_resp.(TTL_name{iO}), 'omitmissing'),length(opto_resp.(TTL_name{iO})))
end
%% plotting Fr and ISI
if plot_flag
    figure

    [~, ~] = MS_kmean_scatter([fr', isi', b_idx'], 2, [1,2,3], 50);
    xlabel('firing rate')
    ylabel('ISI')
    zlabel('bursting index')
    axis square
end

