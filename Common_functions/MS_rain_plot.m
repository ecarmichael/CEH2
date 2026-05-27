function [ax, p, stats_out] =  MS_rain_plot(data_in, group_idx, color, stats, x_vals, marks)

% based off this cool function by Tom Marshall: https://git.fmrib.ox.ac.uk/marshall/public/-/blob/d0b82399a75ea907b332779421c00ff34e2f669a/raincloud_plots/raincloud_plot.m


if nargin < 4
    stats = [];
    x_vals = [1 2]; % where to put the plot. useful for putting multiple together.
    marks = 'off';
elseif nargin < 5
    x_vals = [1 2];
    marks = 'off';
elseif nargin < 6
    marks = 'off';
end
%%
% loop over data groups. 
grps = unique(group_idx); 

for ii = 1:length(grps)
    
hold on

    Y = data_in(group_idx == grps(ii)); 
    
%     cloud dist
[a,b] = ksdensity(Y', 'NumPoints', 200, 'bandwidth', []);

% boxplot
Q = quantile(Y,[0.25 0.75 0.5 0.02 0.98]);


% density plot
% scale_val = (.5 - max(a)) ; 
% a_n = (-a)*scale_val; 
% a_n = (a_n./max(a_n))+.25; 
ax(ii).h{1} = patch(b,((a./max(a))./2)+x_vals(ii)+.05, color(ii,:));

% ax(ii).h{1} = patch(b,(-a)*10*scale_val+x_vals(ii), color(ii,:));

% rotate(ax(ii).h{1}, [0 1 0], 90)
set(ax(ii).h{1}, 'FaceColor', color(ii,:));
set(ax(ii).h{1}, 'EdgeColor', 'none'); %[0.1 0.1 0.1]);
set(ax(ii).h{1}, 'LineWidth', 2);


% box
% data
ax(ii).sc{1} = scatter(Y, (x_vals(ii)-1)+sort(MS_randn_range(length(Y), 1, .65,.95)) ,10,  color(ii,:), 'filled');

% mean line
ax(ii).h{2} = line([Q(3) Q(3)],[0 0.1] + x_vals(ii),'col','k','LineWidth',1);

if strcmpi(marks, 'box')
    ax(ii).box = rectangle('Position',[Q(1) x_vals(ii)-.1 Q(2)-Q(1) .2],'EdgeColor',[0 0 0 .75], 'LineWidth',1.5);

    % ax(ii).h{4} = line([Q(1) Q(1)],[-0.05 0] + x_vals(ii),'col',[0 0 0  .75],'LineWidth',2);
    %
    % ax(ii).h{5} = line([Q(2) Q(2)],[-0.05 0] + x_vals(ii),'col',[0 0 0  .75],'LineWidth',2);

    % whiskers

    ax(ii).eb{1} = line([Q(4) Q(1)], [x_vals(ii) x_vals(ii)],'col',[0 0 0  .5],'LineWidth',1);
    ax(ii).eb{2} = line([Q(2) Q(5)], [x_vals(ii) x_vals(ii)],'col',[0 0 0  .5],'LineWidth',1);
elseif strcmpi(marks, 'wiskers')
    ax(ii).h{6} = line([Q(5) Q(5)], [0 0.05] + x_vals(ii),'col',[0 0 0 .5],'LineWidth',1);
    ax(ii).h{7} = line([Q(4) Q(4)], [0 0.05] + x_vals(ii),'col',[0 0 0  .5],'LineWidth',1);
end


end

view([90 -90])
%%

if ~isempty(stats) 
    if grps > 3
        warning('Skipping STATS: Only set up to do 2 group stats ATM...')
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

    fprintf('| Cohen d = %.2f \n', eff_s.Effect);
    fprintf('Group A (%.2f +/- %.2f)  Vs Group B (%.2f +/- %.2f) \n\n',mean(data_in(group_idx == grps(1)), 'omitnan'),MS_SEM(data_in(group_idx == grps(1))),...
        mean(data_in(group_idx == grps(2)), 'omitnan'),MS_SEM(data_in(group_idx == grps(2))))

end

