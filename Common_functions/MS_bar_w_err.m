function [hb, eb] =  MS_bar_w_err(data_a, data_b, color, xval)

hb = bar([nanmean(data_a), nanmean(data_b)], 'FaceColor', color(1,:), 'EdgeColor', color(1,:));
hold on
eb = errorbar([nanmean(data_a), nanmean(data_b)], [MS_SEM(data_a) ,MS_SEM(data_b)]);
eb.LineStyle = 'none';
eb.Color = 'k';