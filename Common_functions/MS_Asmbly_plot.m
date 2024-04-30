function MS_Asmbly_plot(A_temp, A_cells, map_out, A_loc, fig_dir, info);


if nargin < 5
    fig_dir = [];
    info = []; 
end

%% loop over all the assemblies to visualize

c_ord = MS_linspecer(length(A_cells)); 

figure(3000);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
win = A_loc{1}.win; 


n = ceil(size(A_cells,1)/3);
m = 8;
l = 4;
c_ii = 0;
f_n = 0;
s_plot_max = reshape(1:(l*m), m, l)';

for ii = 1: size(A_cells,2)
    if c_ii+1 > l
        f_n = f_n+1;
        figure(3000+f_n);clf; hold on; set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        c_ii = 1;
    else
        c_ii = c_ii+1;
    end
    
    s_idx = s_plot_max(c_ii,1);
    
    subplot(l,m,[s_idx s_idx+1])
    hold on
    stem(A_temp(:,ii), 'color', [.8 .8 .8 .2])
    stem(A_cells{ii}, A_temp(A_cells{ii},ii), 'color', c_ord(ii,:), 'MarkerFaceColor', c_ord(ii,:))
    view(90,90)
    
    title(['Assembly #' num2str(ii) ' (nPC: ' num2str(sum(map_out{ii}.place_idx))  '/' num2str(length(map_out{ii}.place_idx)) ' = ' num2str(sum(map_out{ii}.place_idx)/length(map_out{ii}.place_idx)*100, '%0.1f') '%)'])
    ylim([-0.1 0.5])
    xlabel('cell ID')
    xlim([0 length(A_temp(:,ii))])
    
    
    subplot(l,m,s_idx+2)
    cla
    hold on
    P_idx = find(logical(map_out{ii}.place_idx));
    
    scatter(A_temp(A_cells{ii},ii),map_out{ii}.MI,35,c_ord(ii,:))
    scatter(A_temp(A_cells{ii}(P_idx),ii),map_out{ii}.MI(P_idx),35,c_ord(ii,:), 'filled')
     ylabel('MI'); xlim([0 0.4]);
    
    % plot the position and trajectory of the animal at the time of the
    % assembly
    subplot(l,m,[s_idx+3])
    cla
    hold on
    for ip = 1:length(A_loc{ii}.loc)
        
        plot(A_loc{ii}.win_time , A_loc{ii}.loc_mat(ip,:), 'color',[c_ord(ii,:) .5], 'linewidth',  3*(A_loc{ii}.peak_val(ip)./max(A_loc{ii}.peak_val)))
    end
    
    
    %         plot((-win:win)/mode(diff(behav.time)), median(this_pos), 'color',[c_ord(ii,:) 1], 'linewidth', 3)
    xlim([A_loc{ii}.win_time(1) A_loc{ii}.win_time(end)]);
    %         set(gca, 'color', 'k')
    plot(0, median(A_loc{ii}.loc), 's','color','k', 'markersize', 20 )
    ylabel('position on track (cm)')
    
    
    subplot(l,m,s_idx+4)
    cla
    hold on
    [N, edges] = histcounts(A_loc{ii}.loc, length(map_out{ii}.bins(1):5:map_out{ii}.bins(end)),'Normalization', 'probability' );
    bar(edges(1:end-1)+mode(diff(edges))/2, N, 'facecolor',[.5 .5 .5], 'faceAlpha', 1, 'EdgeAlpha', 0)

    view(90, 90)
    set(gca,  'xdir', 'reverse'); %'YTick',[0 1], 'yticklabel', {'0' 'max'},
    xlim([map_out{ii}.bins(1) map_out{ii}.bins(end)])
    xlabel('location')
    
    subplot(l,m,s_idx+5)
    cla
    hold on
    [N_C, edges] = histcounts(map_out{ii}.cent*map_out{ii}.bin_size, length(map_out{ii}.bins(1):5:map_out{ii}.bins(end)));
    
    h1 = area(map_out{ii}.bins,map_out{ii}.map_mean./max(map_out{ii}.map_mean,[],2), 'facecolor',c_ord(ii,:), 'EdgeAlpha', 0 );
    h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    h2 = plot(edges(1:end-1)+mode(diff(edges))/2, N_C/max(N_C), 'k', 'linewidth', 1);

    set(gca, 'YTick',[0 1], 'yticklabel', {'0' 'max'}, 'xdir', 'reverse')
    view(90, 90)
   legend('Centroid', 'Location', 'northeast', 'Box', 'off')
    
    
    subplot(l,m,[s_idx+6 s_idx+7])
    
    % sort them based on the centroids
    [~, sort_idx] = sort(map_out{ii}.cent); 
    
    imagesc(1:size(map_out{ii}.map,1),map_out{ii}.bins,  (map_out{ii}.map(sort_idx,:)./max(map_out{ii}.map(sort_idx,:), [],2))')
    hold on
    scatter(1:size(map_out{ii}.cent,1), map_out{ii}.cent(sort_idx)*map_out{ii}.bin_size,55, 'ow');
    scatter(1:size(map_out{ii}.cent,1), map_out{ii}.cent(sort_idx)*map_out{ii}.bin_size,55, 'xw')
    
    set(gca,'xtick',1:size(map_out{ii}.map,1), 'xTickLabel',  A_cells{ii}(sort_idx)','ydir', 'normal', 'XTickLabelRotation', 270);
    ylim([map_out{ii}.bins(1) map_out{ii}.bins(end)])
    xlim([.5 size(map_out{ii}.map,1)+.5])
    colormap(parula(64))
    title(['Z var Peak: ' num2str(round(map_out{ii}.peak_z,2)) '| Cent: ' num2str(round(map_out{ii}.cent_z,2))])
    
    if   c_ii+1 > l
        
        subplot(l,m,[s_idx s_idx+1]);   ylabel('weight')
        
        subplot(l,m,s_idx+2); xlabel('ICA weight');
        
        subplot(l,m,s_idx+3);     xlabel('time from ReAct (s)')
        
        subplot(l,m,s_idx+4);     ylabel('activation(s)')
        
        subplot(l,m,s_idx+5);     ylabel('Mean place map')
        
        subplot(l,m,[s_idx+6 s_idx+7]); xlabel('place cell ID')
        

        if ~isempty(fig_dir)
            saveas(gcf, [fig_dir filesep info.subject '_' info.session '_Wake_A_summery_' num2str(f_n+1) '_' strrep(num2str(info.bin), '.', 'p')  's_bin.png']);
        end
    end
end