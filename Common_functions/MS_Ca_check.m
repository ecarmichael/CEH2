function MS_Ca_check(ms, cell_ids,time, fact)
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
     time = [ms.time(1) ms.time(end)]; 
     
         fact = 1; 
elseif nargin < 3
         time = [ms.time(1) ms.time(end)]; 

        fact = 1; 
elseif nargin < 4
    fact = 1; 
end

if ~isfield(ms, 'time') && isfield(ms, 'tvec')
    ms.time = ms.tvec;
end

%% restrict data 
if (time(1) ~= ms.time(1)) || (time(end) ~= ms.time(end)) 
    fprintf('Restricting between %0.2fs and %0.2fs\n', time(1), time(end))
idx = nearest_idx(time, ms.time ); 

ms.time = ms.time(idx(1):idx(2)); 
ms.deconv = ms.deconv(idx(1):idx(2),:); 
ms.denoise = ms.denoise(idx(1):idx(2),:); 
ms.RawTraces = ms.RawTraces(idx(1):idx(2),:); 
ms.Binary = ms.Binary(idx(1):idx(2),:); 


end
%%
figure(1919)
    clf

% c_ord = [linspace(0,1,length(cell_ids))', zeros(length(cell_ids),2)]
% colormap(colorMap);
    c_ord = parula(length(cell_ids)*2);
    
    subplot(2,3,3)
    cla
    if isfield(ms, 'SFPs_sharp')
    MS_plot_all_SFPs(ms.SFPs_sharp); 
%     imagesc(max(ms.SFPs_sharp(:,:,:),[],3));
%     c_val = caxis; 
%     caxis([0 c_val(2)*.2])

    hold on

    for ii = 1:length(cell_ids)
        [row, col] = find(ms.SFPs_sharp(:,:,cell_ids(ii)) == max(max(ms.SFPs_sharp(:,:,cell_ids(ii)))));
        text(col, row, num2str(cell_ids(ii)), 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'g');
%         scatter(col, row,50,'o',  'MarkerEdgeColor',c_ord(ii,:) , 'LineWidth', 1);% c_ord(ii,:)
    end
    title(['nCells: ' num2str(ms.numNeurons)]);
%     xlim([min(ms.Centroids(:,1))-min(ms.Centroids(:,1))*.2  max(ms.Centroids(:,1))+max(ms.Centroids(:,1))*.2])
%     ylim([min(ms.Centroids(:,2))-min(ms.Centroids(:,2))*.2  max(ms.Centroids(:,2))+max(ms.Centroids(:,2))*.2])
    end
    subplot(2,3,6)
    cla
        binsize = 0.1; % in seconds, so everything else should be seconds too
        gauss_window = 1./binsize; % 1 second window
        gauss_SD = 0.5./binsize; % 0.02 seconds (20ms) SD
        gk = gausskernel(gauss_window,gauss_SD); gk = gk./binsize; % normalize by binsize
        gau_sdf = conv2((ms.deconv./ms.denoise)>0.01,gk,'same'); % convolve with gaussian window
        gau_z = zscore(gau_sdf, [], 2);
        
        
    imagesc( ms.time,1:ms.numNeurons,  gau_z);  set(gca, 'YDir', 'normal'); 
      c_val = caxis; 
    caxis([0 c_val(2)*.2])
%     caxis([0 .01])

%     all_cord = parula(floor(ms.numNeurons*1.2));%zeros(length(ms.units),3);%repmat(linspecer(5), 100, 1); 
%     hold on
%     mult_fac = .05; 
%     for ii = 1:ms.numNeurons
%         this_spk = find((ms.deconv(:,ii)./ms.denoise(:,ii))>0.01);
%         plot([ms.time(this_spk) ms.time(this_spk)]', [ones(length(this_spk),1)*ii-0.5 ones(length(this_spk),1)*ii+0.5]', 'color', all_cord(ii,:))
% %         plot(ms.time, zscore(ms.denoise(:,ii))+ii*mult_fac,  'color', all_cord(ii,:), 'linewidth', 1.5)
% 
% %         plot(ms.frame*(1/30), zscore(ms.Deconv(:,ii))+ii*10,'color', all_cord(ii,:), 'linewidth', 1)
% %         plot(ms.frame*(1/30), (ms.Spikes(:,ii)*5)+ii*10,'color', 'k')
%         
%     end
%        set(gca, 'color', 'k')
    title('zCsp traces')

    ylabel('cell ID')
    xlabel('time (s)')
     c_val = caxis; 
    caxis([0 16])
%     set(gca,'ytick', 0:100:ms.numNeurons*mult_fac, 'YTickLabel', (0:100:length(ms.units)*mult_fac)/mult_fac, 'TickDir', 'out')
    ylim([0 ms.numNeurons])
%     xlim([ms.time(1) ms.time(end)])
%         set(gca,'xtick', [round(abs(ms.time(1))) round(ms.time(end),0)], 'xTickLabel', [round(abs(ms.time(1))) round(ms.time(end),0)], 'TickDir', 'out')

    
    
    subplot(2,3,[1 2 4 5])
    cla
    hold on
    for ii = 1:length(cell_ids)
        this_deconv = ms.deconv(:,cell_ids(ii)); 
        this_deconv(this_deconv==0) = NaN; 

        plot(ms.time, (ms.RawTraces(:,cell_ids(ii)))+ii*fact,'color', [0.6 0.6 0.6], 'linewidth', 1)
        plot(ms.time, (ms.denoise(:,cell_ids(ii)))+ii*fact,'color', [c_ord(ii,:) .5], 'linewidth', 1)
        plot(ms.time, ((this_deconv*5) - fact/4)+ii*fact,'color','k', 'linewidth',1)
        if isfield(ms, 'Binary')
            this_binary = ms.Binary(:,cell_ids(ii))/10; 
             this_binary(this_binary>1) = NaN; 
        plot(ms.time, ((this_binary)- fact/4)+ii*fact,'color','k', 'linewidth',1)
        end
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
    xlim([round(abs(ms.time(1))) round(ms.time(end),0)])
    legend({'Raw', 'Denoise', 'Deconv'}, 'Location', 'north', 'Orientation', 'horizontal', 'box', 'off')
    

% maximize
