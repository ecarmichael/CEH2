function MS_asmbly_quick_plot(temp, proj,data_in,  idx)
%% MS_asmbly_quick_plot: plots the templates and projection for a specific assembly given the template and projection arrays




subplot(2,4,[1 5])

hold on
stem(temp(:,idx), 'color', [.8 .8 .8 .2])

a_idx = sum(zscore(temp(:,idx)) > 2, 2) > 0;

stem(find(a_idx), temp(find(a_idx),idx), 'color',winter(1), 'MarkerFaceColor', winter(1))

view(90,90)

ax(1) = subplot(2,4,2:4);

plot(proj(idx, :));


ax(2) = subplot(2,4,6:8);

imagesc(data_in(:, a_idx)')

linkaxes(ax, 'x')

end