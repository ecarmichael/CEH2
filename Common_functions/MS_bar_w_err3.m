function [hb, eb, sc, p, stats_out] =  MS_bar_w_err3(data_a, data_b,data_c, color, data_flag, stats, x_vals)


if nargin < 4
    data_flag = 1; 
    stats = 'anova2'; 
    x_vals = [1 2 3]; % where to put the plot. useful for putting multiple together.
    color = MS_linspecer(3); 
elseif nargin < 5
    data_flag = 1; 
    stats = 'anova2'; 
    x_vals = [1 2 3]; % where to put the plot. useful for putting multiple together. 
elseif nargin < 6
    stats = 'anova2'; 
    x_vals = [1 2 3];
elseif nargin < 7
    x_vals = [1 2 3];
end

if data_flag == 1
    hb = bar(x_vals, [mean(data_a, 'omitnan'), mean(data_b, 'omitnan'), mean(data_c, 'omitnan')]', 'FaceColor', 'flat', 'EdgeColor','flat');
    hb.CData(1,:) = color(1,:); 
    hb.CData(2,:) = color(2,:); 
    hb.CData(3,:) = color(3,:); 

end
hold on

offsets_a = x_vals(1)+ sort(MS_randn_range(length(data_a), 1, -.1, .1));
offsets_b = x_vals(2)+ sort(MS_randn_range(length(data_b), 1, -.1, .1));
offsets_c = x_vals(3)+ sort(MS_randn_range(length(data_c), 1, -.1, .1));




if ~isempty(stats)
    % zero-pad for uneven Ns
d_mat = NaN(3, max([size(data_a,2),size(data_b,2), size(data_c,2)]));
d_mat(1,1:length(data_a)) = data_a; 
d_mat(2,1:length(data_b)) = data_b; 
d_mat(3,1:length(data_c)) = data_c;
% make a table
tbl = table(d_mat(1,:)', d_mat(2,:)', d_mat(3,:)', 'Variablenames', {'A', 'B', 'C'}); 
meas = table([1 2 3]', 'VariableNames',{'meas'});

    switch stats
        case 'anova1'
            % disp('using ttest')
            [p,stats_out.a_tbl, stats_out.stats] = anova1([d_mat(1,:)', d_mat(2,:)', d_mat(3,:)'],[], 'off');
            stats_out.stats.F = stats_out.a_tbl{2,5}; 

            if p(1) > 0.05; h =1; else h=0; end

            [results, means] = multcompare(stats_out.stats, 'display', 'off');
            stats_out.m_tbl =  array2table([results,means],"VariableNames", ...
                ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value","Mean","Standard Error"]);

        case 'anova2'
            % disp('using ttest')
            [p, stats_out.a_tbl, stats_out.stats] = anova2([d_mat(1,:)', d_mat(2,:)', d_mat(3,:)']', 3, "off");
            stats_out.stats.F = stats_out.a_tbl{2,5}; 
            
            if p(1) < 0.05; h =1; else h=0; end

            [results, means] = multcompare(stats_out.stats,  'display', 'off');
            stats_out.m_tbl =  array2table([results,means],"VariableNames", ...
                ["Group A","Group B","Lower Limit","A-B","Upper Limit","P-value","Mean","Standard Error"]);

        plot([offsets_a, offsets_b]', [data_a ;data_b;], '-', 'Color', [.5 .5 .5])
        plot([offsets_b, offsets_c]', [data_b ;data_c;], '-', 'Color', [.5 .5 .5])


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

    % print the results
    if p(1) < 0.05
        fprintf('<strong>%s</strong> - F(<strong>%d</strong>,<strong>%d</strong>): <strong>%.2f</strong> p = <strong>%.5f</strong>\n',stats, stats_out.a_tbl{2,3},stats_out.a_tbl{3,3}, stats_out.a_tbl{2,5}, stats_out.a_tbl{2,6})
    else
        fprintf('<strong>%s</strong> - F(%d,%d): %.2f p = %.5f\n',stats, stats_out.a_tbl{2,3},stats_out.a_tbl{3,3}, stats_out.a_tbl{2,5})
    end
    % overall
    if ~isnan(h) && p(1) < 0.05
        if size(data_a,2) == 1
            data_pool = [data_a; data_b; data_c]; 
        else
            data_pool = [data_a, data_b, data_c];
        end
        ylim([0 max(data_pool)*1.3])
        x_lim = xlim; 
        if (0.05 > p(1)) && (p(1)> 0.01)
            text(x_lim(end)*.9, max(data_pool, [], 'all')*1.3, ['*overall p = ' num2str(p(1), 3)], 'color', 'k', 'FontSize',10, 'HorizontalAlignment','right')
        elseif (0.01 >= p(1)) && (p(1) >= 0.001)
            text(x_lim(end), max(data_pool, [], 'all')*1.3, ['** overall p = ' num2str(p(1), 3)], 'color', 'k', 'FontSize',10, 'HorizontalAlignment','right')
        elseif p(1) < 0.001
            text(x_lim(end), max(data_pool, [], 'all')*1.3, ['*** overall p = ' num2str(p(1), 3)], 'color', 'k', 'FontSize',10, 'HorizontalAlignment','right')
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
            text(median(x_vals(1:2)), max(data_pool, [], 'all')*1.15, ['* p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',10)
        elseif (0.01 >= this_p) && (this_p >= 0.001)
            text(median(x_vals(1:2))*.975, max(data_pool, [], 'all')*1.15, ['** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',10)
        elseif this_p < 0.001
            text(median(x_vals(1:2))*.975, max(data_pool, [], 'all')*1.15, ['*** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',10)
        end

        plot(x_vals(1:2), [max(data_pool, [], 'all')*1.1 max(data_pool, [], 'all')*1.1], '-k', 'linewidth', 1)
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
            text(median([x_vals(1) x_vals(3)]), max(data_pool, [], 'all')*1.05, ['* p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',10)
        elseif (0.1 >= this_p) && (this_p >= 0.001)
            text(median([x_vals(1) x_vals(3)])*.975, max(data_pool, [], 'all')*1.05, ['** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',10)
        elseif this_p < 0.001
            text(median([x_vals(1) x_vals(3)])*.95, max(data_pool, [], 'all')*1.05, ['*** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',10)
        end

        plot([x_vals(1) x_vals(3)], [max(data_pool, [], 'all')*1.025 max(data_pool, [], 'all')*1.025], '-k', 'linewidth', 1)
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
            text(median(x_vals(2:3)), max(data_pool, [], 'all')*1.1, ['* p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',10)
        elseif (0.1 >= this_p) && (this_p >= 0.001)
            text(median(x_vals(2:3))*.975, max(data_pool, [], 'all')*1.1, ['** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',10)
        elseif this_p < 0.001
            text(median(x_vals(2:3))*.95, max(data_pool, [], 'all')*1.1, ['*** p = ' num2str(this_p, 3)], 'color', 'k', 'FontSize',10)
        end

        plot(x_vals(2:3), [max(data_pool, [], 'all')*1.05 max(data_pool, [], 'all')*1.05], '-k', 'linewidth', 1)
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

end

% add in the data points (after the lines) so they are the top layer
if data_flag > 0 %&& length(data_a) == length(data_b)
        sc{1} = scatter( offsets_a, data_a,15, 'markerfacecolor', 'k', 'MarkerEdgeColor', 'none');
    sc{2} = scatter(offsets_b, data_b,15,  'markerfacecolor', 'k', 'MarkerEdgeColor', 'none');
        sc{3} = scatter(offsets_c, data_c,15,  'markerfacecolor', 'k', 'MarkerEdgeColor', 'none');

    % sc{1} = scatter( offsets_a, data_a,25, 'markerfacecolor', color(1,:), 'MarkerEdgeColor', [.2 .2 .2]);
    % sc{2} = scatter(offsets_b, data_b,25,  'markerfacecolor', color(2,:), 'MarkerEdgeColor', [.2 .2 .2]);
else
    sc = []; 
end

eb = errorbar(x_vals, [mean(data_a, 'omitnan'), mean(data_b,'omitnan'), mean(data_c,'omitnan')], [MS_SEM(data_a) ,MS_SEM(data_b), MS_SEM(data_c)]);
eb.LineStyle = 'none';
eb.Color = [.2 .2 .2];
eb.LineWidth =1; 

if ~exist('hb','var')
    hb = sc;
end

