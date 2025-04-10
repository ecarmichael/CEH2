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
    hb = bar(x_vals, [nanmean(data_a), nanmean(data_b)], 'FaceColor', color(1,:), 'EdgeColor', color(1,:));
end
hold on


if data_flag > 0 %&& length(data_a) == length(data_b)
    sc{1} = scatter(1+sort(MS_randn_range(length(data_a), 1, -.1, .1)), data_a,25,  color(1,:), 'filled');
    sc{2} = scatter(2+sort(MS_randn_range(length(data_b), 1, -.1, .1)), data_b,25,  color(2,:), 'filled');
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
        case 'ttest2'
            disp('using ttest2')
            [h, p, ~,stats] = ttest2(data_a, data_b);
            
        case 'ranksum'
            disp('using ranksum')
            [h, p, ~,stats] = ranksum(data_a, data_b);
    end
    
    if h
        plot(median(x_vals), max([data_a, data_b], [], 'all')*1.15, '*', 'color', 'k')
        plot(x_vals, [max([data_a, data_b], [], 'all')*1.1 max([data_a, data_b], [], 'all')*1.1], '-k', 'linewidth', 1.5)
    end
end

if ~exist('hb','var')
    hb = sc;
end

