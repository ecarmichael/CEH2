function MS_Asmbly_plot_raster_figure(A_in, fig_dir,plot_idx, type)

if nargin < 2
    fig_dir = [];
    plot_idx = 1:length(A_in.P_pos);
    type = 'act_mat';
elseif nargin < 3
    plot_idx = 1:length(A_in.P_pos);
    type = 'act_mat';
elseif nargin < 4
    type = 'act_mat';
end

%%
ft_size = 21
% c_ord = MS_linspecer(length(plot_idx)+ceil(length(plot_idx)/5));
c_ord = [67, 127, 151; 132 147 35; 255 179 13; 253 22 26]/255;
% make a colour coded activity array
act_array = zeros(size(A_in.wake_data));

cnt = 0; c_list = []; c_map = [0 0 0]; c_step = 0;
for iA = plot_idx
    c_idx = A_in.P_pos{iA};
    A_range = 1:length(c_idx);
    A_range = A_range+cnt;
    cnt = cnt+length(A_range)+1;
    
    for ii = 1:length(c_idx)
        act_array(logical(A_in.wake_data(:,c_idx(ii))),A_range(ii)) = find(iA ==plot_idx);
    end
    
    
    c_list = [c_list; c_idx];
    c_map = [c_map; c_ord(find(iA ==plot_idx),:)];
    c_step = [c_step cnt];
    c_label{find(iA ==plot_idx)} = [num2str(iA)];
end


%%


f = figure(109); 
% f.WindowState = 'maximized';
clf
set(gcf, 'units','inches', 'position',[0 0 6 4], 'innerposition',[0 0 6 4],'outerposition',[0 0 6 4])
% set(gcf,  'PaperPosition', [0 0 6 3])
% axes('Units', 'normalized', 'Position', [0 0 1 1])

% maximize
pause(.5)

ax(1) = subplot(3, 1, 1);
cla
hold on
plot(A_in.behav.time/1000, A_in.behav.position(:,1), '-', 'color', [.5 .5 .5], 'linewidth', 0.3)
scatter(A_in.behav.time(A_in.move_idx)/1000, A_in.behav.position((A_in.move_idx),1),ones(size(A_in.behav.time(A_in.move_idx)))*1, A_in.behav.speed(A_in.move_idx), 'Marker', '.')
%     plot(A_in.behav.time(A_in.move_idx)/1000, A_in.behav.position((A_in.move_idx),1));
xlim([min(A_in.behav.time(A_in.move_idx)) max(A_in.behav.time(A_in.move_idx))]/1000)
% ylabel({'position (cm)'})
set(gca, 'XTick', [], 'linewidth', 1, 'ytick', [0 100]);

if isstr(type) == 1
    ax(2) = subplot(3,1,2);
    cla
    hold on
    off_set = 0; these_idx = [];
    
    if strcmpi(type, 'act_mat')
        imagesc(A_in.wake_tvec, 1:size(act_array,2),  act_array')
        ylim([0.5 length(c_list)+0.5])
%         set(gca, 'YTick', c_step(2:end) - diff(c_step)/2, 'YTickLabel',c_label, 'xticklabel', [])
        set(gca, 'XTickLabel', [], 'yticklabels', []); 
        h = get(gca);
        h.XAxis.TickLength = [ 0 0];
        
        for ii = 2:length(c_step)
           text(0, c_step(ii) - diff(c_step(ii-1:ii))/2, c_label{ii-1},'HorizontalAlignment', 'right', 'Color', c_map(ii,:), 'FontSize', ft_size )
        end
        
    else
        
        for iA = plot_idx
            
            for ii = size(A_in.P_pos{iA}, 1):-1:1
                iC = A_in.P_pos{iA}(ii);
                off_set = off_set+1;
                
                if strcmpi(type, 'raw')
                    plot(A_in.wake_tvec, ((A_in.wake_data(:,iC)./max(A_in.wake_data(:,iC)))*.8)+off_set,  'color', c_ord(iA,:), 'linewidth', 1)
                elseif strcmpi(type, 'raster')
                    this_idx = find(A_in.wake_data(:,iC) > 0);
                    this_t = A_in.wake_tvec(this_idx);
                    this_s = (ones(size(A_in.wake_tvec(this_idx)))*ii);
                    plot([this_t; this_t], [this_s-.5+off_set; this_s+.5+off_set], 'color', c_ord(iA,:), 'linewidth', 2)
                    %               plot(rec.time, rec.Binary(:,iC)*ii +.5+off_set, 'color', c_ord(iA,:), 'linewidth', 2)
                    %             plot([S.t{iC}, S.t{iC}]', [(ones(size(S.t{iC}))*ii)-.5+off_set, (ones(size(S.t{iC}))*ii)+.5+off_set]', 'color', c_ord(iA,:), 'linewidth', 2)
                end
            end
            these_idx = [these_idx, A_in.P_pos{iA}']; % keep track to avoid overlap;
        end
        ylim([0.5 off_set+0.5])
        
    end
    
else
    ax(2) = subplot(7,1,3:3);
    cla
    hold on
    
    imagesc(A_in.wake_tvec, 1:size(act_array,2),  act_array')
    ylim([0 length(c_list)])
    set(gca, 'YTick', c_step(2:end) - diff(c_step)/2, 'YTickLabel',c_label)
    
    
    ax(4) = subplot(7,1,4:5);
    cla
    hold on
    off_set = 0; these_idx = [];
    for iA = plot_idx
        
        for ii = size(A_in.P_pos{iA}, 1):-1:1
            iC = A_in.P_pos{iA}(ii);
            off_set = off_set+1;
            
            plot(A_in.wake_tvec, ((A_in.wake_data(:,iC)./max(A_in.wake_data(:,iC)))*.8)+off_set,  'color', c_ord(find(iA == plot_idx),:), 'linewidth', .5)
            
        end
        these_idx = [these_idx, A_in.P_pos{iA}']; % keep track to avoid overlap;
    end
    ylim([0 length(c_list)+1])
    
end
%%
%     non_ass_idx = 1:length(S.t);
%     rm_idx = (ismember(non_ass_idx, these_idx));
%     non_ass_idx(rm_idx) = [];
%
%     for ii = size(non_ass_idx, 2):-1:1
%         if ~isempty(S.t{non_ass_idx(ii)})
%             plot([S.t{non_ass_idx(ii)}, S.t{non_ass_idx(ii)}]', [(ones(size(S.t{non_ass_idx(ii)}))*ii)-.5+off_set, (ones(size(S.t{non_ass_idx(ii)}))*ii)+.5+off_set]', 'color', [.7 .7 .7 .7], 'linewidth', 2)
%         end
%     end
%     ylim([0 size(S.t, 2)])
%     set(gca, 'XTick', [], 'YDir', 'normal', 'ytick', [])
%     ylabel('Cell activity')

ax(3) = subplot(3,1,3);
cla
hold on
leg_val = [];
yline(log10(10), '--', 'color', [.7 .7 .7], 'linewidth', 0.3)

for iA = plot_idx
    this_proj = A_in.P_proj(iA,:);
    %         this_proj(this_proj <=0) =NaN;
    plot(A_in.wake_tvec, log10(this_proj), 'color', c_ord(find(iA == plot_idx),:), 'linewidth', .5)
    leg_val{find(iA == plot_idx)} = [ num2str(iA)];
end
ylim([0 inf])
ylabel({'assembly strength'})
xlabel('time (s)')
y_val = get(gca, 'YTick'); 
set(gca, 'YTickLabel', 10.^y_val, 'linewidth', 1); 

legend(['R thresh'], 'Orientation', 'horizontal', 'box', 'off')



linkprop(ax,{'XLim'});
colormap(ax(1), 'parula');
colormap(ax(2),[1 1 1; c_map(2:end,:)]);
if length(ax) == 4
    colormap(ax(4),c_map);
    
end

%% save the figure
cfg.ft_size= 21;
cfg.resize = 0;
SetFigure(cfg, gcf)
% set(gcf,'PaperOrientation','landscape');
% 
% set(gcf,'PaperUnits','inches', 'Units', 'inches');
% set(gcf,'PaperSize', [1.5 1]);

if ~isempty(fig_dir)
    %         saveas(gcf, [fig_dir filesep A_in.info.subject '_' A_in.info.session '_' strrep(num2str(A_in.info.bin), '.', 'p') 's_bin_wake_' type '.png']);
%     if exist([fig_dir filesep A_in.info.subject '_' A_in.info.session '_' strrep(num2str(A_in.info.bin), '.', 'p') 's_bin_wake_raster.pdf'], 'file')
%         
%         delete([fig_dir filesep A_in.info.subject '_' A_in.info.session '_' strrep(num2str(A_in.info.bin), '.', 'p') 's_bin_wake_raster.pdf'], 'file');
%     end
    print(gcf, '-dpdf', [fig_dir filesep A_in.info.subject '_' A_in.info.session '_' strrep(num2str(A_in.info.bin), '.', 'p') 's_bin_wake_raster.pdf'])
    linkprop(ax,{'XLim'});

    xlim([200 300])
    ax(2) = subplot(3,1,2);
    for ii = 2:length(c_step)
        text(200-4, c_step(ii) - diff(c_step(ii-1:ii))/2, c_label{ii-1},'HorizontalAlignment', 'right', 'Color', c_map(ii,:), 'FontSize', ft_size )
    end
        print(gcf, '-dpdf', [fig_dir filesep A_in.info.subject '_' A_in.info.session '_' strrep(num2str(A_in.info.bin), '.', 'p') 's_bin_wake_raster_zoom.pdf'])

end

%% format for figure panels




