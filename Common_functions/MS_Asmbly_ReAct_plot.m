function MS_Asmbly_ReAct_plot(A_out, fig_dir, type)

if nargin < 3
    type = 'REM';
end

c_ord = MS_linspecer(size(A_out.P_proj,1));


x_val = [min([A_out.([type '_Pre_tvec'])(1)  A_out.([type '_Post_tvec'])(1)]), max([A_out.([type '_Pre_tvec'])(end)  A_out.([type '_Post_tvec'])(end)])];

A_idx = 1:size(A_out.P_temp,2);
figure(310);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

ar(1)= subplot(10,1,1:2);
imagesc(A_out.([type '_Pre_tvec']),1:size(A_out.([type '_Pre_data']),2),   A_out.([type '_Pre_data'])')
%     xlim([tvec_rem_pre(1) tvec_rem_pre(end)])
xlim(x_val)
title([A_out.info.subject ' ' A_out.info.session])
set(gca, 'xtick', [])


ar(2) = subplot(10,1,3:4);
cla;
hold on
for ii = A_idx
    if A_out.([type '_Pre_stats']).p_val(ii) <0.05
        plot(A_out.([type '_Pre_tvec']), A_out.([type '_Pre_proj'])(A_idx(ii),:), 'color', c_ord(ii,:), 'linewidth', 1)
    else
        plot(A_out.([type '_Pre_tvec']), A_out.([type '_Pre_proj'])(A_idx(ii),:),'--', 'color', c_ord(ii,:),'linewidth', 0.5)
    end
end
ylim([0 50])
%     xlim([tvec_rem_pre(1) tvec_rem_pre(end)])
xlim(x_val)
ylabel(['Pre ' type ' Reactivation'])
yline(A_out.([type '_Pre_stats']).R_thresh)
legend(num2str(A_idx'),'location',  'northeast', 'box', 'off', 'orientation', 'vertical', 'NumColumns', ceil(length(A_idx)/2))


ap(1) = subplot(10,1,5:6);
imagesc(A_out.([type '_Post_tvec']),1:size(A_out.([ type '_Post_data']),2),   A_out.([type '_Post_data']))
%     xlim([tvec_rem(1) tvec_rem(end)])
xlim(x_val)
set(gca, 'xtick', [])

ap(2) = subplot(10,1,7:8);
cla;
hold on
for ii = A_idx
    if A_out.([ type '_Post_stats']).p_val(ii) <0.05
        plot(A_out.([ type '_Post_tvec']), A_out.([ type '_Post_proj'])(A_idx(ii),:), 'color', c_ord(ii,:), 'linewidth', 1)
    else
        plot(A_out.([ type '_Post_tvec']), A_out.([ type '_Post_proj'])(A_idx(ii),:),'--', 'color', c_ord(ii,:),'linewidth', 0.5)
    end
end
ylim([0 50])
%     xlim([tvec_rem(1) tvec_rem(end)])
xlim(x_val)
ylabel(['Post ' type ' Reactivation'])
yline(A_out.([ type '_Post_stats']).R_thresh)
%     legend(num2str(A_idx'),'location',  'northeast', 'box', 'off', 'orientation', 'horizontal')

if size(A_out.([ type '_Post_shuff']).data,2) > size(A_out.([ type '_Pre_shuff']).data,2)
    s_phase = 'Post';
else
    s_phase = 'Pre';
end

as(1) = subplot(10,1,9);
c = imagesc(A_out.([ type '_' s_phase '_tvec']),1:size(A_out.([ type '_' s_phase '_shuff']).data,2),   A_out.([ type '_' s_phase '_shuff']).data');

xlim([A_out.([ type '_' s_phase '_tvec'])(1) A_out.([ type '_' s_phase '_tvec'])(end)])
set(gca, 'xtick', [])

g_ord = [linspace(0.3, .9, length(A_idx)) ; linspace(0.3, .9, length(A_idx)) ; linspace(0.3, .9, length(A_idx))]';
as(2) = subplot(10,1,10);
cla;
hold on
for ii = A_idx
    plot(A_out.([ type '_' s_phase '_tvec']), A_out.([ type '_' s_phase '_shuff']).proj(ii,:), 'color', g_ord(find(ii == A_idx),:))
end
ylim([0 50])
%     xlim([tvec_rem(1) tvec_rem(end)])
xlim(x_val)
ylabel(['Shuff ' s_phase ' Reactivation'])
xlabel('time (s)')
yline(A_out.([ type '_' s_phase '_stats']).R_thresh)

linkaxes(ar, 'x');
linkaxes(ap, 'x');
linkaxes(as, 'x');

if fig_dir
    saveas(gcf, [fig_dir filesep A_out.info.subject '_' A_out.info.session '_' strrep(num2str(A_out.info.bin), '.', 'p') 's_bin_' type '_A_summary.png']);
end

end