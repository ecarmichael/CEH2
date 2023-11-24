function h = SWR_plot_CA_nlx(ms, csc, cell_ids)
%% SWR_plot_CA_nlx:
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
% EC 2023-11-23   initial version 
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
figure(1901)
    clf


    c_ord = parula(length(cell_ids)*2);
    
 ax(1)= subplot(6,1,1:5)
    cla
    hold on
    fact = 2.5; 
    for ii = 1:length(cell_ids)
        this_deconv = ms.deconv(:,cell_ids(ii)); 
        this_deconv(this_deconv==0) = NaN; 
        plot(ms.time, (ms.RawTraces(:,cell_ids(ii)))+ii*fact,'color', [0.6 0.6 0.6], 'linewidth', 1)
        plot(ms.time, (ms.denoise(:,cell_ids(ii)))+ii*fact,'color', [c_ord(ii,:) .5], 'linewidth', 1)
        plot(ms.time, ((this_deconv*5) - fact/4)+ii*fact,'color','k', 'linewidth',1)
        if isfield(ms, 'Binary')
        plot(ms.time, ((this_deconv*5) - fact/4)+ii*fact,'color','k', 'linewidth',1)
        end
%         plot(ms.time, ((this_deconv*5) - fact/4)+ii*fact,'color', c_ord(ii,:), 'linewidth', 1)
%         plot(ms.time, (ms.Spikes(:,ii)*5)+ii*10,'color', 'k')
        
    end
    title('Denoise and deconvolved traces for sample cells')
    ylabel('cell ID')
%     xlabel('time (s)')
    set(gca,'ytick', [1 length(cell_ids)]*fact, 'YTickLabel', [cell_ids(1) cell_ids(end)], 'TickDir', 'out')
%     set(gca, 'xticklabel', []); 

    ylim([0 (length(cell_ids)+2)*fact])
    xlim([round(abs(ms.time(1))) round(ms.time(end),0)])
    legend({'Raw', 'Denoise', 'Deconv'}, 'Location', 'north', 'Orientation', 'horizontal', 'box', 'off')
    
   ax(2) =  subplot(6,1,6);
    plot(csc.tvec - csc.tvec(1), csc.data(1,:)); 
    xlim([0 csc.tvec(end) - csc.tvec(1)])
%             set(gca,'xtick', [round(abs(ms.time(1))) round(ms.time(end),0)], 'xTickLabel', [round(abs(ms.time(1))) round(ms.time(end),0)], 'TickDir', 'out')

linkaxes(ax, 'x')    

% maximize


