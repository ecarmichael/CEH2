function [cfs,frq, tvec_data, freq]= MS_plot_cw_emg(tvec, csc, emg, idx,  win_s, cent_tvec)
%% MS_plot_cw_emg: plot the continuous wavelet transform along with the EMG signal centered on 'idx' index. 
%
%
%
%    Inputs: 
%    - tvec: [1 x nSamples]  time vector for LFP and EMG
%    - csc: [1 x nSamples]   LFP data for the channel of interest.
%
%    - emg: [1 x nSamples]   EMG data corresponding to the LFP data. 
%
%    - win_s: [2 x 1 double]   how many seconds before and after the event
%
%    - cent_tvec: [logical]   swith timing to be centered on event. 
%
%    Outputs: 
%    - h:  [figure handle]
%
%
%
%
% EC 2022-05-24   initial version 
%
%
%
%% initialize
if nargin < 6
    cent_tvec = 0
end

Fs = round(1/mode(diff(tvec))); 

%%
        CSC_data = csc(1,idx - (win_s*Fs):idx + (win_s*Fs));
        
        if ~isempty(emg) % prep EMG if you have it.
            EMG_data = emg(1,idx - (win_s*Fs):idx + (win_s*Fs));
        end
        cwt(CSC_data,  Fs);
        [cfs,frq] = cwt(CSC_data,  Fs);
        xlabels = get(gca, 'xtick'); 
        set(gca, 'xticklabels', xlabels - win_s)
        h = gcf; 
        x_lim = xlim;
        hold on
        xline(win_s, '--w', 'start', 'linewidth', 2);
%         xline(x_lim(2) - win_s, '--k', 'end', 'linewidth', 2);
        yline(10, '--w', '10hz', 'linewidth', 2, 'LabelHorizontalAlignment','left');
        
        AX = gca;
        [minf,maxf] = cwtfreqbounds(numel(CSC_data),Fs);
        
        freq = 2.^(round(log2(minf)):round(log2(maxf)));
        AX.YTickLabelMode = 'auto';
        AX.YTick = freq;
        ylim([1 140])
        
        %add the LFP
        tvec_data = tvec(idx - (win_s*Fs):idx + (win_s*Fs),1);
        plot(tvec_data - tvec_data(1), (CSC_data*1600)+4, 'w')
        
        % add EMG if you have it.
        if ~isempty(emg)
            plot(tvec_data - tvec_data(1), (EMG_data*1200)+2, 'color', [.7 .7 .7])
        end
        
%         pause(2); % 
        
%%add in raster plot.  
%         figure(iR*100);
%         % restrict the spike times to this event.
%         S_r = restrict(S, pREM_IV.tstart(iR)-win_s, pREM_IV.tend(iR)+win_s);
%         for iS = 1:length(S_r.t)
%             S_r.t{iS} = S_r.t{iS} - (pREM_IV.tstart(iR)-win_s); % offset for plotting.
%         end
%         
%         % make a raster.
%         cfg_ras = [];
%         cfg_ras.LineWidth = 1;
%         cfg_ras.spkColor = 'linspecer';
%         cfg_ras.openNewFig = 0;
%         MultiRaster(cfg_ras, S_r)
%         %        c_mat = [linspecer(sum(place_idx));  repmat([1 1 1],sum(~place_idx),1)]; % make colors depending on thecentroid
%         xline(win_s, '--w', 'start', 'linewidth', 2);
%         x_lim = xlim;
%         xline(x_lim(2) - win_s, '--w', 'end', 'linewidth', 2);
%         set(gca, 'color', 'k'); %set background color.
%         %         colormap([linspecer(sum(place_idx));  repmat([1 1 1],1,1)]); % used for centroids
%         %         cx = colorbar;
%         %         cx.TickLabels = cx.Ticks * max(centroids);
%         %         cx.Label.String = 'place cell centroid';
%         
%         
% %         pause(1)
%         % move both plots to new plot.
%         figlist=get(groot,'Children');
%         
%         newfig=figure(iR*1000);
%         tcl=tiledlayout(2,1);
%         
%         for jj= 1:numel(figlist)
%             if (figlist(jj).Number ~= iR*10) && (figlist(jj).Number ~= iR*100)
%                 continue
%             end
%             figure(figlist(jj));
%             ax=gca;
%             ax.Parent=tcl;
%             ax.Layout.Tile=jj;
%         end
%         close(iR*10);
%         close(iR*100);
%         set(gcf, 'position', [662 96 758 892])
%         set(gcf, 'InvertHardcopy', 'off')
%         
%         pause(1)