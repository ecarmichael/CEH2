function [hb, eb, sc, p, stats] =  MS_bar_w_err(data_a, data_b, color, data_flag, stats, x_vals)

if nargin < 4
    data_flag = 0; 
    stats = []; 
    x_vals = [1 2]; % where to put the plot. useful for putting multiple together. 
elseif nargin < 5
    stats = []; 
    x_vals = [1 2];
elseif nargin < 6
    x_vals = [1 2];
end

if data_flag == 1
    hb = bar(x_vals, [nanmean(data_a), nanmean(data_b)]', 'FaceColor', 'flat');
    hb.CData(1,:) = color(1,:); 
    hb.CData(2,:) = color(2,:); 
end
hold on

offsets_a = x_vals(1)+ sort(MS_randn_range(length(data_a), 1, -.1, .1));
offsets_b = x_vals(2)+ sort(MS_randn_range(length(data_b), 1, -.1, .1));

if data_flag > 0 %&& length(data_a) == length(data_b)
    sc{1} = scatter( offsets_a, data_a,25, 'markerfacecolor', color(1,:), 'MarkerEdgeColor', [.2 .2 .2]);
    sc{2} = scatter(offsets_b, data_b,25,  'markerfacecolor', color(2,:), 'MarkerEdgeColor', [.2 .2 .2]);
else
    sc = []; 
end

eb = errorbar(x_vals, [nanmean(data_a), nanmean(data_b)], [MS_SEM(data_a) ,MS_SEM(data_b)]);
eb.LineStyle = 'none';
eb.Color = 'k';
eb.LineWidth =1.5; 

if ~isempty(stats)
    switch stats
        case 'ttest'
            disp('using ttest')
            [h, p, ~, stats] = ttest(data_a, data_b);

            % add connections between points. 

            plot([offsets_a, offsets_b]', [data_a ;data_b], '-', 'Color', [.6 .6 .6])
        case 'ttest2'
            disp('using ttest2')
            [h, p, ~,stats] = ttest2(data_a, data_b);
            
        case 'ranksum'
            disp('using ranksum')
            [h, p, ~,stats] = ranksum(data_a, data_b);
    end
    
    if ~isnan(h) && h == 1
        if size(data_a,2) == 1
            data_pool = [data_a; data_b]; 
        else
            data_pool = [data_a, data_b];
        end
        if (0.5 > p) && (p > 0.01)
            text(median(x_vals), max(data_pool, [], 'all')*1.1, '*', 'color', 'k', 'FontSize',22)
        elseif (0.1 >= p) && (p >= 0.001)
            text(median(x_vals)*.975, max(data_pool, [], 'all')*1.1, '**', 'color', 'k', 'FontSize',22)
        elseif p < 0.001
            text(median(x_vals)*.95, max(data_pool, [], 'all')*1.1, '***', 'color', 'k', 'FontSize',22)
        end
                   
        text(median(x_vals)*1.1, max(data_pool, [], 'all')*1.15, ['   p = ' num2str(p, 3)], 'color', 'k', 'FontSize',16)

        plot(x_vals, [max(data_pool, [], 'all')*1.05 max(data_pool, [], 'all')*1.05], '-k', 'linewidth', 1.5)
    end
end

if ~exist('hb','var')
    hb = sc;
end

