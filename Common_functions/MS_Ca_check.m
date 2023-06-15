function MS_Ca_check(ms, cell_ids)
%% MS_Ca_check: Quick plot to check traces, deconv, and SFPs in ms file. 
%
%
%
%    Inputs: 
%    -
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2023-03-21   initial version 
%
%
%
%% initialize  

if nargin < 2
    cell_ids = 1:20; 
     if ms.numNeurons < 20
        cell_ids  = 1:ms.numNeurons;
     end
end

if ~isfield(ms, 'time') && isfield(ms, 'tvec')
    ms.time = ms.tvec;
end

%%
figure(1919)
    clf

    c_ord = summer(length(cell_ids));
    
    subplot(2,2,2)
    cla
    MS_plot_all_SFPs(ms.SFPs_sharp); 
%     imagesc(max(ms.SFPs_sharp(:,:,:),[],3));
%     c_val = caxis; 
%     caxis([0 c_val(2)*.2])

    hold on

    for ii = 1:length(cell_ids)
        [row, col] = find(ms.SFPs_sharp(:,:,cell_ids(ii)) == max(max(ms.SFPs_sharp(:,:,cell_ids(ii)))));
        scatter(col, row,50,'o',  'MarkerEdgeColor',c_ord(ii,:) , 'LineWidth', 1);% c_ord(ii,:)
    end
    title(['nCells: ' num2str(ms.numNeurons)]);
%     xlim([min(ms.Centroids(:,1))-min(ms.Centroids(:,1))*.2  max(ms.Centroids(:,1))+max(ms.Centroids(:,1))*.2])
%     ylim([min(ms.Centroids(:,2))-min(ms.Centroids(:,2))*.2  max(ms.Centroids(:,2))+max(ms.Centroids(:,2))*.2])

    subplot(2,2,4)
    cla
    imagesc( ms.time,1:ms.numNeurons,  zscore(ms.Binary)'); set(gca, 'YDir', 'normal'); 
      c_val = caxis; 
    caxis([0 c_val(2)*.2])
%     all_cord = parula(floor(ms.numNeurons*1.2));%zeros(length(ms.units),3);%repmat(linspecer(5), 100, 1); 
%     hold on
%     mult_fac = .05; 
%     for ii = 1:ms.numNeurons
%         plot(ms.time, zscore(ms.denoise(:,ii))+ii*mult_fac,  'color', all_cord(ii,:), 'linewidth', 1.5)
% 
% %         plot(ms.frame*(1/30), zscore(ms.Deconv(:,ii))+ii*10,'color', all_cord(ii,:), 'linewidth', 1)
% %         plot(ms.frame*(1/30), (ms.Spikes(:,ii)*5)+ii*10,'color', 'k')
%         
%     end
%        set(gca, 'color', 'k')
    title('zscored deconvolved traces')

    ylabel('cell ID')
    xlabel('time (s)')
     c_val = caxis; 
    caxis([0 16])
%     set(gca,'ytick', 0:100:ms.numNeurons*mult_fac, 'YTickLabel', (0:100:length(ms.units)*mult_fac)/mult_fac, 'TickDir', 'out')
    ylim([0 ms.numNeurons])
%     xlim([ms.time(1) ms.time(end)])
xlim([65 260])
%         set(gca,'xtick', [round(abs(ms.time(1))) round(ms.time(end),0)], 'xTickLabel', [round(abs(ms.time(1))) round(ms.time(end),0)], 'TickDir', 'out')

    
    
    subplot(2,2,[1 3])
    cla
    hold on
    fact = 2.5; 
    for ii = 1:length(cell_ids)
        this_deconv = ms.deconv(:,cell_ids(ii)); 
        this_deconv(this_deconv==0) = NaN; 
        plot(ms.time, (ms.RawTraces(:,cell_ids(ii)))+ii*fact,'color', [0.6 0.6 0.6], 'linewidth', 1)
        plot(ms.time, (ms.denoise(:,cell_ids(ii)))+ii*fact,'color', [c_ord(ii,:) .5], 'linewidth', 1)
        plot(ms.time, ((this_deconv*5) - fact/4)+ii*fact,'color','k', 'linewidth',1)

%         plot(ms.time, ((this_deconv*5) - fact/4)+ii*fact,'color', c_ord(ii,:), 'linewidth', 1)
%         plot(ms.time, (ms.Spikes(:,ii)*5)+ii*10,'color', 'k')
        
    end
    title('Denoise and deconvolved traces for sample cells')
    ylabel('cell ID')
    xlabel('time (s)')
    set(gca,'ytick', [1 length(cell_ids)]*fact, 'YTickLabel', [cell_ids(1) cell_ids(end)], 'TickDir', 'out')
%     x_ticks = get(gca, 'xtick'); 
%         x_lab = get(gca, 'xticklabel'); 
        set(gca,'xtick', [round(abs(ms.time(1))) round(ms.time(end),0)], 'xTickLabel', [round(abs(ms.time(1))) round(ms.time(end),0)], 'TickDir', 'out')

    ylim([0 (length(cell_ids)+2)*fact])
    xlim([ms.time(1) ms.time(end)])
    legend({'Raw', 'Denoise', 'Deconv'}, 'Location', 'north', 'Orientation', 'horizontal', 'box', 'off')
    

maximize
