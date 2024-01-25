function MS_Asmbly_ReAct_plot(A_out, fig_dir); 

c_ord = MS_linspecer(size(A_out.P_proj,1)); 


x_val = [min([A_out.REM_Pre_tvec(1)  A_out.REM_Post_tvec(1)]), max([A_out.REM_Pre_tvec(end)  A_out.REM_Post_tvec(end)])];

    A_idx = 1:size(A_out.P_temp,2);
    figure(310);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    ar(1)= subplot(10,1,1:2);
    imagesc(A_out.REM_Pre_tvec,1:size(A_out.REM_Pre_data,2),   A_out.REM_Pre_data')
    %     xlim([tvec_rem_pre(1) tvec_rem_pre(end)])
    xlim(x_val)
    title([A_out.info.subject ' ' A_out.info.session])
    set(gca, 'xtick', [])
    
    
    ar(2) = subplot(10,1,3:4);
    cla;
    hold on
    for ii = A_idx
        if A_out.REM_Pre_stats.p_val(ii) <0.05
            plot(A_out.REM_Pre_tvec, A_out.REM_Pre_proj(A_idx(ii),:), 'color', c_ord(ii,:), 'linewidth', 1)
        else
            plot(A_out.REM_Pre_tvec, A_out.REM_Pre_proj(A_idx(ii),:),'--', 'color', c_ord(ii,:),'linewidth', 0.5)
        end
    end
    ylim([0 50])
    %     xlim([tvec_rem_pre(1) tvec_rem_pre(end)])
    xlim(x_val)
    ylabel('Pre REM Reactivation')
    yline(A_out.REM_Pre_stats.R_thresh)
    legend(num2str(A_idx'),'location',  'northeast', 'box', 'off', 'orientation', 'vertical', 'NumColumns', ceil(length(A_idx)/2))
    
    
    ap(1) = subplot(10,1,5:6);
    imagesc(A_out.REM_Post_tvec,1:size(A_out.REM_Post_data,2),   A_out.REM_Post_data')
    %     xlim([tvec_rem(1) tvec_rem(end)])
    xlim(x_val)
    set(gca, 'xtick', [])
    
    ap(2) = subplot(10,1,7:8);
    cla;
    hold on
    for ii = A_idx
        if A_out.REM_Post_stats.p_val(ii) <0.05
            plot(A_out.REM_Post_tvec, A_out.REM_Post_proj(A_idx(ii),:), 'color', c_ord(ii,:), 'linewidth', 1)
        else
            plot(A_out.REM_Post_tvec, A_out.REM_Post_proj(A_idx(ii),:),'--', 'color', c_ord(ii,:),'linewidth', 0.5)
        end
    end
    ylim([0 50])
    %     xlim([tvec_rem(1) tvec_rem(end)])
    xlim(x_val)
    ylabel('Post REM Reactivation')
    yline(A_out.REM_Post_stats.R_thresh)
%     legend(num2str(A_idx'),'location',  'northeast', 'box', 'off', 'orientation', 'horizontal')
    
    
    as(1) = subplot(10,1,9);
    c = imagesc(A_out.REM_Post_tvec,1:size(A_out.REM_Post_shuff.data,2),   A_out.REM_Post_shuff.data');
    
    xlim([A_out.REM_Post_tvec(1) A_out.REM_Post_tvec(end)])
    set(gca, 'xtick', [])
    
    g_ord = [linspace(0.3, .9, length(A_idx)) ; linspace(0.3, .9, length(A_idx)) ; linspace(0.3, .9, length(A_idx))]';
    as(2) = subplot(10,1,10);
    cla;
    hold on
    for ii = A_idx
        plot(A_out.REM_Post_tvec, A_out.REM_Post_shuff.proj(ii,:), 'color', g_ord(find(ii == A_idx),:))
    end
    ylim([0 50])
    %     xlim([tvec_rem(1) tvec_rem(end)])
    xlim(x_val)
    ylabel('Shuff postReactivation')
    xlabel('time (s)')
    yline(A_out.REM_Post_stats.R_thresh)
    
    linkaxes(ar, 'x');
    linkaxes(ap, 'x');
    linkaxes(as, 'x');
    
    if fig_dir
        saveas(gcf, [fig_dir filesep A_out.info.subject '_' A_out.info.session '_' strrep(num2str(A_out.info.bin), '.', 'p') 's_bin_REM_A_summary.png']);
    end
    
end