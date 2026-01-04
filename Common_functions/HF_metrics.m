function [data_out] = HF_metrics(data_in, dat_path)
%% HF_Metrics:  computes simple spiking measures

if nargin < 2
    dat_path = [];
end


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
    c_red = [0.9153    0.2816    0.2878];
    TTL = {'6', '7'};
    TTL_name = {'red', 'blue'};

    for iO = 1:2
        if sum(contains(data_in.evts.label, TTL{iO}))
            this_TTL_idx = find(contains(data_in.evts.label, TTL{iO}));

            evt_t = data_in.evts.t{this_TTL_idx} ;

            % isolate events of a certain length
            e_d = evt_t(2,:) - evt_t(1,:);
            ITIs = unique(round(e_d, 3));

            for iTi = length(ITIs):-1:1

                evt_t = data_in.evts.t{this_TTL_idx} ;
                k_idx = round(e_d,3) == ITIs(iTi) ;
                evt_t(:,~k_idx) = [];

                % isolate the cell
                S = data_in.S;
                S.t = [];
                S.t{1} = data_in.S.t{iS};
                S.label = [];
                S.label{1} =  data_in.S.label{iS};

                evt_t = sort(evt_t(1,:));

                cfg_peth = [];
                cfg_peth.window = [-.25 .25];
                cfg_peth. plot_type = 'zscore';
                cfg_peth.dt = 0.001;
                cfg_peth.gauss_window = .025;
                cfg_peth.gauss_sd = .0025;
                cfg_peth.shuff = 500;
                cfg_peth.t_on = mode(e_d(k_idx));
                cfg_peth.rec_color = c_red;
                cfg_peth.plot = 'on';
                [peth_S, ~,~, ~, ~, ~, ~, ~, ~,peth_T] = SpikePETH_Shuff(cfg_peth, S, evt_t);

                % pre_fr = mean(zval(peth_IT<0));
                % laser_fr = mean(zval(find(peth_IT==0):nearest_idx(mode(e_d(k_idx)), peth_IT)));

                % t-test for pre and post rate for each stim?
                u_val = unique(peth_T); pre = NaN(1,length(u_val)); post = pre;
                for iV = length(u_val):-1:1
                    this_idx = peth_T == u_val(iV);
                    if isempty(this_idx)
                        continue
                    end
                    pre(iV) = sum((peth_S(this_idx) > -mode(e_d(k_idx))) & (peth_S(this_idx) < 0));
                    post(iV) = sum((peth_S(this_idx) > 0) & peth_S(this_idx) < mode(e_d(k_idx)));
                end

                S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.dur = mode(e_d(k_idx));
                S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.ITI = mode(diff(evt_t));
                S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.pre = pre;
                S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.post = post;
                % if ~isempty(pre) || ~isempty(post)
                [~,S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.pval] = ttest(pre, post);

                if iTi == 1
                    opto_resp.(TTL_name{iO})(iS) = ttest(pre, post);
                end
                % else
                % S_metrics{iS}.(['opto_' TTL_name{iO}]){iTi}.pval = NaN;
                % opto_resp.(TTL_name{iO})(iS) = NaN;
            end
        end % end ITIS
    end
end% opto colours

end

%%
for ii = 1:length(S_metrics)
    
    fprintf('Red Response Cell #%.d, %.4f\n',ii,  S_metrics{ii}.opto_red{1}.pval)

end

fprintf('Red optoresponsive cells %.2f%% (%.d/%.d)\n', (nansum(opto_resp.red)./length(opto_resp.red))*100, nansum(opto_resp.red),length(opto_resp.red))
fprintf('Red optoresponsive cells %.2f%% (%.d/%.d)\n', (nansum(opto_resp.red)./length(opto_resp.red))*100, nansum(opto_resp.red),length(opto_resp.red))

%% plotting Fr and ISI

figure(9990)

[g_idx, n_val] = MS_kmean_scatter([fr', isi', b_idx'], 2, [1,2,3], 50);
xlabel('firing rate')
ylabel('ISI')
zlabel('bursting index')
axis square
