function [hb, eb, p] =  MS_bar_w_err(data_a, data_b, color, data_flag, stats, x_vals)

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

hb = bar(x_vals, [nanmean(data_a), nanmean(data_b)], 'FaceColor', color(1,:), 'EdgeColor', color(1,:));
hold on


if data_flag && length(data_a) == length(data_b)
    for ii = length(data_a):-1:1
        plot(x_vals(1), data_a(ii), '.', 'color', [.5 .5 .5], 'markersize', 15)
        plot(x_vals(1), data_b(ii), '.', 'color', [.5 .5 .5], 'markersize', 15)
        plot(x_vals, [data_a(ii) data_b(ii)], '-', 'color', [.5 .5 .5], 'linewidth', .5)
    end
end

eb = errorbar(x_vals, [nanmean(data_a), nanmean(data_b)], [MS_SEM(data_a) ,MS_SEM(data_b)]);
eb.LineStyle = 'none';
eb.Color = 'k';
eb.LineWidth =1.5; 

if ~isempty(stats)
    switch stats
        case 'ttest'
            disp('using ttest')
            [h, p] = ttest(data_a, data_b);
        case 'ttest2'
            disp('using ttest2')
            [h, p] = ttest2(data_a, data_b);
            
        case 'ranksum'
            disp('using ranksum')
            [h, p] = ranksum(data_a, data_b);
    end
    
    if h
        plot(median(x_vals), max([data_a, data_b], [], 'all')*1.15, '*', 'color', 'k')
        plot(x_vals, [max([data_a, data_b], [], 'all')*1.1 max([data_a, data_b], [], 'all')*1.1], '-k', 'linewidth', 1.5)
    end
end

