function [ax, p, stats_out] =  MS_bar_n(data_in, group_idx, color, stats, x_vals, marks)

% based off this cool function by Tom Marshall: https://git.fmrib.ox.ac.uk/marshall/public/-/blob/d0b82399a75ea907b332779421c00ff34e2f669a/raincloud_plots/raincloud_plot.m

error('WIP')

if nargin < 4
    stats = [];
    x_vals = []; % where to put the plot. useful for putting multiple together.
    marks = 'off';
elseif nargin < 5
    x_vals = [];
    marks = 'off';
elseif nargin < 6
    marks = 'off';
end

% loop over data groups.
grps = unique(group_idx);

if isempty(color)
g_cord = MS_linspecer(length(grps)); % if colour is not specified. 
else
    g_ord  = color; 
end

if isempty(x_vals)
    x_vals = 1:numel(grps);
else
    x_vals = x_vals(:)';
end

if isempty(color)
    color = MS_linspecer(numel(grps));
end
%% loop over groups. 
mean_vals = NaN(size(grps)); 
sem_vals = NaN(size(grps));
n_vals = sem_vals; 

cla
hold on
for ii = 1:length(grps)
    mean_vals(ii) = mean(data_in(group_idx == grps(ii)), 'omitnan'); 
    sem_vals(ii) = MS_SEM(data_in(group_idx == grps(ii))); 
    n_vals(ii) = sum(group_idx == grps(ii)); 


    hb = bar(x_vals(ii), mean_vals(ii),'FaceColor', 'flat', 'EdgeColor','flat');
    hb.CData(ii,:) = color(ii,:);

    offsets = x_vals(ii)+ sort(MS_randn_range(length(data_in(group_idx == grps(ii))), 1, -.1, .1)); 
    sc{ii} = scatter( offsets, data_in(group_idx == grps(ii)),15, 'markerfacecolor', 'k', 'MarkerEdgeColor', 'none');

end


% errorbars
eb = errorbar(x_vals, mean_vals, sem_vals);
eb.LineStyle = 'none';
eb.Color = [.2 .2 .2 .2];
eb.LineWidth =1;
eb(1).CapSize = eb(1).CapSize*2;

if ~exist('hb','var')
    hb = sc;
end

%%  make a data array 
        d_mat = NaN(length(grps), max(n_vals));
        var_names = []; 
        for ii = 1:length(grps)
            d_mat(ii,1:n_vals(ii)) = data_in(group_idx == grps(ii));
            var_names{ii} = num2str(ii); 
        end

%%

    if length(grps) == 3

        switch stats
            case 'anova1'
                [p,stats_out.a_tbl, stats_out.stats] = anova1(d_mat,[], 'off');
                stats_out.stats.F = stats_out.a_tbl{2,5};


            if p(1) > 0.05; h =1; else h=0; end

            [results, means] = multcompare(stats_out.stats, 'display', 'off');
            stats_out.stats.std = means(:,2);

            stats_out.m_tbl =  array2table(results,"VariableNames", ...
                ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value"]);


            case 'anova2'
                [p, stats_out.a_tbl, stats_out.stats] = anova2(d_mat, size(d_mat), "off");
                stats_out.stats.F = stats_out.a_tbl{2,5};

                if p(1) < 0.05; h =1; else, h=0; end

                [results, means] = multcompare(stats_out.stats,  'display', 'off');
                stats_out.m_tbl =  array2table([results,means],"VariableNames", ...
                    ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value","Mean","Standard Error"]);

                plot([offsets_a, offsets_b]', [data_a ,data_b]', '-', 'Color', [.5 .5 .5])
                plot([offsets_b, offsets_c]', [data_b ,data_c]', '-', 'Color', [.5 .5 .5])


                % case 'ranova'
                %     % disp('using ttest2')
                %     rm = fitrm(tbl, 'A-C ~ 1', 'WithinDesign',meas);
                %     [stats_out.r_tbl, ~, ~, h] = ranova(rm);
                %
                %     % reconstruct to match ANOVA
                %
                %     results = multcompare(rm,'meas');
                %     % add connections between points.
                %     plot([offsets_a, offsets_b, offsets_c]', [data_a ;data_b; data_c], '-', 'Color', [.5 .5 .5])
            case 'KW' % nonparametric version of anova1
                % disp('using ks')
                [h, p, stats_out.stats] = kruskalwallis([d_mat(1,:)', d_mat(2,:)', d_mat(3,:)']);
        end

        % 3 way ploting
    if p(1) < 0.05
        fprintf('<strong>%s</strong> - F(<strong>%d</strong>,<strong>%d</strong>): <strong>%.2f</strong> p = <strong>%.5f</strong>\n',stats, stats_out.a_tbl{2,3},stats_out.stats.df, stats_out.stats.F, stats_out.a_tbl{2,6})
    else
        fprintf('<strong>%s</strong> - F(%d,%d): %.2f p = %.5f \n',stats, stats_out.a_tbl{2,3},stats_out.stats.df, stats_out.stats.F, stats_out.a_tbl{2,6})
    end
    % overall
    if ~isnan(h) && p(1) < 0.05
        if size(data_a,2) == 1
            data_pool = [data_a; data_b; data_c];
        else
            data_pool = [data_a, data_b, data_c];
        end
        xlim([0 max(data_pool)*1.3])
        y_lim = ylim;
        % add the sig markers
        if (0.05 > p(1)) && (p(1)> 0.01)
            text(max(data_pool, [], 'all')*1.25, y_lim(end)*.9, ['*overall p = ' num2str(p(1), 3)], 'color', 'k', 'FontSize',8, 'HorizontalAlignment','right')
        elseif (0.01 >= p(1)) && (p(1) >= 0.001)
            text( max(data_pool, [], 'all')*1.25,y_lim(end),y_lim(end), ['** overall p = ' num2str(p(1), 3)], 'color', 'k', 'FontSize',8, 'HorizontalAlignment','right')
        elseif p(1) < 0.001
            text( max(data_pool, [], 'all')*1.25,y_lim(end), ['*** overall p = ' num2str(p(1), 3)], 'color', 'k', 'FontSize',8, 'HorizontalAlignment','right')
        elseif p(1) < 0.0001
            text( max(data_pool, y_lim(end),[], 'all')*1.25, ['**** overall p = ' num2str(p(1), 3)], 'color', 'k', 'FontSize',8, 'HorizontalAlignment','right')
        end
    end

    % 1 vs 2
    this_p = stats_out.m_tbl.("P-value")(1);

    if ~isnan(this_p) && this_p < 0.05
        if size(data_a,2) == 1
            data_pool = [data_a; data_b];
        else
            data_pool = [data_a, data_b];
        end
        if (0.05 > this_p) && (this_p> 0.01)
            text( max(data_pool, [], 'all')*1.15, median(x_vals(1:2)),['* p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',8)
        elseif (0.01 >= this_p) && (this_p >= 0.001)
            text( max(data_pool, [], 'all')*1.15,median(x_vals(1:2))*.975, ['** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',8)
        elseif this_p < 0.001
            text( max(data_pool, [], 'all')*1.15,median(x_vals(1:2))*.975, ['*** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',8)
        elseif this_p < 0.0001
            text( max(data_pool, [], 'all')*1.15,median(x_vals(1:2))*.975, ['**** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',8)
        end

        plot( [max(data_pool, [], 'all')*1.1 max(data_pool, [], 'all')*1.1],x_vals(1:2), '-k', 'linewidth', 1)
    end

    % 1 vs 3
    this_p = stats_out.m_tbl.("P-value")(2);

    if ~isnan(this_p) && this_p < 0.05
        if size(data_a,2) == 1
            data_pool = [data_a; data_c];
        else
            data_pool = [data_a, data_c];
        end
        if (0.5 > this_p) && (this_p> 0.01)
            text( max(data_pool, [], 'all')*1.05,median([x_vals(1) x_vals(3)]), ['* p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',8)
        elseif (0.1 >= this_p) && (this_p >= 0.001)
            text( max(data_pool, [], 'all')*1.05,median([x_vals(1) x_vals(3)])*.975, ['** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',8)
        elseif this_p < 0.001
            text( max(data_pool, [], 'all')*1.05, median([x_vals(1) x_vals(3)])*.95,['*** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',8)
        elseif this_p < 0.0001
            text( max(data_pool, [], 'all')*1.05, median([x_vals(1) x_vals(3)])*.925,['**** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',8)
        end

        plot( [max(data_pool, [], 'all')*1.025 max(data_pool, [], 'all')*1.025], [x_vals(1) x_vals(3)],'-k', 'linewidth', 1)
    end


    % 2 vs 3
    this_p = stats_out.m_tbl.("P-value")(3);

    if ~isnan(this_p) && this_p < 0.05
        if size(data_a,2) == 1
            data_pool = [data_b; data_c];
        else
            data_pool = [data_b, data_c];
        end
        if (0.5 > this_p) && (this_p> 0.01)
            text( max(data_pool, [], 'all')*1.1, median(x_vals(2:3)),['* p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',8)
        elseif (0.1 >= this_p) && (this_p >= 0.001)
            text( max(data_pool, [], 'all')*1.1, median(x_vals(2:3))*.975,['** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',8)
        elseif this_p < 0.001
            text(max(data_pool, [], 'all')*1.1,median(x_vals(2:3))*.95,  ['*** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',8)
        elseif this_p < 0.0001
            text( max(data_pool, [], 'all')*1.1,median(x_vals(2:3))*.95, ['**** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',8)
        end

        plot([max(data_pool, [], 'all')*1.05 max(data_pool, [], 'all')*1.05],x_vals(2:3),  '-k', 'linewidth', 1)
    end

    % report the posthoc
    if isfield(stats_out, 'm_tbl')
        if stats_out.m_tbl{1,6} < 0.05
            fprintf('Group A (%.2f +/- %.2f)  Vs Group B (%.2f +/- %.2f); p =  <strong>%.5f</strong> \n',mean(data_a, 'omitnan'),MS_SEM(data_a),mean(data_b, 'omitnan'),MS_SEM(data_b),  stats_out.m_tbl{1,6})
        else
            fprintf('Group A (%.2f +/- %.2f)  Vs Group B (%.2f +/- %.2f); p =  %.5f \n',mean(data_a, 'omitnan'),MS_SEM(data_a),mean(data_b, 'omitnan'),MS_SEM(data_b),  stats_out.m_tbl{1,6})
        end
        if stats_out.m_tbl{2,6} < 0.05
            fprintf('Group A (%.2f +/- %.2f)  Vs Group C (%.2f +/- %.2f); p =  <strong>%.5f</strong> \n',mean(data_a, 'omitnan'),MS_SEM(data_a),mean(data_c, 'omitnan'),MS_SEM(data_c),  stats_out.m_tbl{2,6})
        else
            fprintf('Group A (%.2f +/- %.2f)  Vs Group C (%.2f +/- %.2f); p =  %.5f \n',mean(data_a, 'omitnan'),MS_SEM(data_a),mean(data_c, 'omitnan'),MS_SEM(data_c),  stats_out.m_tbl{2,6})
        end
        if stats_out.m_tbl{3,6} < 0.05
            fprintf('Group B (%.2f +/- %.2f)  Vs Group C (%.2f +/- %.2f); p =  <strong>%.5f</strong> \n',mean(data_b, 'omitnan'),MS_SEM(data_b),mean(data_c, 'omitnan'),MS_SEM(data_c),  stats_out.m_tbl{3,6})
        else
            fprintf('Group B (%.2f +/- %.2f)  Vs Group C (%.2f +/- %.2f); p =  %.5f \n',mean(data_b, 'omitnan'),MS_SEM(data_b),mean(data_c, 'omitnan'),MS_SEM(data_c),  stats_out.m_tbl{3,6})
        end
    end

    fprintf('\n'); %extra spacing.

    else
        switch stats
            case 'ttest'
                disp('using ttest')
                [h, p,~, stats_out] = ttest(data_in(group_idx == grps(1)), data_in(group_idx == grps(2)));
                eff_s = meanEffectSize(data_in(group_idx == grps(1)), data_in(group_idx == grps(2)), "Effect","cohen", 'Paired',true);

            case 'ttest2'
                disp('using ttest2')
                [h, p,~, stats_out] = ttest2(data_in(group_idx == grps(1)), data_in(group_idx == grps(2)));
                eff_s = meanEffectSize(data_in(group_idx == grps(1)), data_in(group_idx == grps(2)), "Effect","cohen");

            case 'ranksum'
                disp('using ranksum')
                [h, p,~, stats_out] = ranksum(data_in(group_idx == grps(1)), data_in(group_idx == grps(2)));
                eff_s = meanEffectSize(data_in(group_idx == grps(1)), data_in(group_idx == grps(2)), "Effect","cohen");
        end
        if h
            x_lim = xlim;
            if (0.5 > p) && (p > 0.01)
                text(max(data_in)+((x_lim(2)-x_lim(1))*.25),median(x_vals), '*', 'color', 'k', 'FontSize',10)
            elseif (0.1 >= p) && (p >= 0.001)
                text(max(data_in)+((x_lim(2)-x_lim(1))*.25),median(x_vals)*.95, '**', 'color', 'k', 'FontSize',10)
            elseif p < 0.001
                text(max(data_in)+((x_lim(2)-x_lim(1))*.25),median(x_vals)*.925, '***', 'color', 'k', 'FontSize',10)
            elseif p < 0.0001
                text(max(data_in)+((x_lim(2)-x_lim(1))*.25),median(x_vals)*.9, '****', 'color', 'k', 'FontSize',10)
            end
            % plot(max(data_in)+((x_lim(2)-x_lim(1))*.25),median(x_vals),  '*', 'color', 'k')
            plot([max(data_in)+((x_lim(2)-x_lim(1))*.2) max(data_in)+((x_lim(2)-x_lim(1))*.2)], x_vals, '-k', 'linewidth', 1.5)
        end

    if p(1) < 0.05
        fprintf('<strong>%s</strong> - t(<strong>%d</strong>) = <strong>%.2f</strong>, p = <strong>%.4f</strong> ',stats, stats_out.df,stats_out.tstat, p)
    else
        fprintf('<strong>%s</strong> - t(%d) = %.2f, p = %.5f ',stats, stats_out.df,stats_out.tstat, p)
    end

    fprintf('| Cohen d = %.3f \n', eff_s.Effect);
    fprintf('Group A (%.3f +/- %.3f)  Vs Group B (%.3f +/- %.3f) \n\n',mean(data_in(group_idx == grps(1)), 'omitnan'),MS_SEM(data_in(group_idx == grps(1))),...
        mean(data_in(group_idx == grps(2)), 'omitnan'),MS_SEM(data_in(group_idx == grps(2))))

    end


end

