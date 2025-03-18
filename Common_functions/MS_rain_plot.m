function [ax, p] =  MS_rain_plot(data_in, group_idx, color, stats, x_vals)

% based off this cool function by Tom Marshall: https://git.fmrib.ox.ac.uk/marshall/public/-/blob/d0b82399a75ea907b332779421c00ff34e2f669a/raincloud_plots/raincloud_plot.m


if nargin < 4
        stats = []; 
    x_vals = [1 2]; % where to put the plot. useful for putting multiple together. 
elseif nargin < 5
    x_vals = [1 2];

end
%%
% loop over data groups. 
grps = unique(group_idx); 
cla
for ii = 1:size(grps, 1)
    
hold on

    Y = data_in(group_idx == grps(ii)); 
    
%     cloud dist
[a,b] = ksdensity(Y', 'NumPoints', 200, 'bandwidth', []);

% boxplot
Q = quantile(Y,[0.25 0.75 0.5 0.02 0.98]);


% density plot
scale_val = (.5 - max(a)) ; 
ax(ii).h{1} = patch(b,(-a)*10*scale_val+x_vals(ii), color(ii,:));

% rotate(ax(ii).h{1}, [0 1 0], 90)
set(ax(ii).h{1}, 'FaceColor', color(ii,:));
set(ax(ii).h{1}, 'EdgeColor', color(ii,:)); %[0.1 0.1 0.1]);
set(ax(ii).h{1}, 'LineWidth', 2);


% box

% mean line
ax(ii).h{3} = line([Q(3) Q(3)],[0.1 0] + x_vals(ii),'col','k','LineWidth',2);
% h{2} = rectangle('Position',[Q(1) ii Q(2)-Q(1) .5]);
ax(ii).h{4} = line([Q(1) Q(1)],[0.05 0] + x_vals(ii),'col',[color(ii,:) 1],'LineWidth',2);

ax(ii).h{5} = line([Q(2) Q(2)],[0.05 0] + x_vals(ii),'col',[color(ii,:) 1],'LineWidth',2);

% whiskers
ax(ii).h{6} = line([Q(5) Q(5)], [0.05 0] + x_vals(ii),'col',[color(ii,:) .5],'LineWidth',2);
ax(ii).h{7} = line([Q(4) Q(4)],[0.05 0] + x_vals(ii),'col',[color(ii,:) .5],'LineWidth',2);
    
% data
ax(ii).sc{1} = scatter(Y, x_vals(ii)+sort(MS_randn_range(length(Y), 1, 0,.2)) +.1 ,25,  color(ii,:), 'filled');
    

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
                [h, p] = ttest(data_in(group_idx == grps(1)), data_in(group_idx == grps(2)));
            case 'ttest2'
                disp('using ttest2')
                [h, p] = ttest2(data_in(group_idx == grps(1)), data_in(group_idx == grps(2)));
                
            case 'ranksum'
                disp('using ranksum')
                [h, p] = ranksum(data_in(group_idx == grps(1)), data_in(group_idx == grps(2)));
        end
    end
    
    if h
        x_lim = xlim; 
        
        plot(max(data_in)+((x_lim(2)-x_lim(1))*.25),median(x_vals),  '*', 'color', 'k')
        plot([max(data_in)+((x_lim(2)-x_lim(1))*.2) max(data_in)+((x_lim(2)-x_lim(1))*.2)], x_vals, '-k', 'linewidth', 1.5)
    end
end

