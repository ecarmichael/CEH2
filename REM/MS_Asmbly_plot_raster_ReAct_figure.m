function MS_Asmbly_plot_raster_ReAct(A_in, fig_dir,sleep_phase,plot_idx, c_ord , type)


    
    if strcmpi(sleep_phase, 'pREM_data')
        REM_tvec = [ 'REM_Pre_tvec'];
        REM_proj = ['pREM_proj']; 
        sleep_phase = 'REM_Pre_data'; 
        A_in.P_pos = A_in.pREM_A_pos;
    else
            REM = strsplit(sleep_phase, 'data'); 
    REM_tvec = [ REM{1} 'tvec'];
    REM_proj = [REM{1} 'proj']; 
    end
    



if nargin < 3
    fig_dir = [];
    plot_idx = 1:length(A_in.P_pos);
    type = 'act_mat';
            c_ord = []; 

elseif nargin < 4
    plot_idx = 1:length(A_in.P_pos);
    type = 'act_mat';
            c_ord = []; 

    elseif nargin < 5
    type = 'act_mat';
        c_ord = []; 
elseif nargin < 6
    type = 'act_mat';
end



%%

ft_size = 14; 
if isempty(c_ord);
    if length(plot_idx) <= 4
        c_ord = [67, 127, 151; 132 147 35; 255 179 13; 253 22 26]/255;
    else
        c_ord = MS_linspecer(length(plot_idx)+ceil(length(plot_idx)/5));
    end
end
    
   
% make a colour coded activity array
act_array = zeros(size(A_in.(sleep_phase)));

cnt = 0; c_list = []; c_map = [1 1 1]; c_step = 0;
for iA = plot_idx
    c_idx = A_in.P_pos{iA};
    A_range = 1:length(c_idx);
    A_range = A_range+cnt;
    cnt = cnt+length(A_range);
    
    for ii = 1:length(c_idx)
        act_array(logical(A_in.(sleep_phase)(:,c_idx(ii))),A_range(ii)) = find(iA ==plot_idx);
    end
    
    
    c_list = [c_list; c_idx];
    c_map = [c_map; c_ord(find(iA ==plot_idx),:)];
    c_step = [c_step cnt];
    c_label{find(iA ==plot_idx)} = [num2str(iA)];
end



    %%
    f = figure(110); 
% f.WindowState = 'maximized';
clf
set(gcf, 'units','inches', 'position',[0 0 1.5 1]*6, 'innerposition',[0 0 1.5 1]*6,'outerposition',[0 0 1.5 1]*6)
% set(gcf, 'units','inches', 'position',[0 0 1.5 1]*6, 'innerposition',[0 0 1.5 1]*6,'outerposition',[0 0 1.5 1]*6)
% set(gcf,  'PaperPosition', [0 0 6 3])
% axes('Units', 'normalized', 'Position', [0 0 1 1])

% maximize


if isstr(type) == 1
ax(1) = subplot(3, 1, 2);
    cla
    hold on
    off_set = 0; these_idx = [];
    
    if strcmpi(type, 'act_mat')
        imagesc(A_in.(REM_tvec), 1:size(act_array,2),  act_array')
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
                    plot(A_in.(REM_tvec), ((A_in.(sleep_phase)(:,iC)./max(A_in.(sleep_phase)(:,iC)))*.8)+off_set,  'color', c_ord(iA,:), 'linewidth', 1)
                elseif strcmpi(type, 'raster')
                    this_idx = find(A_in.(sleep_phase)(:,iC) > 0);
                    this_t = A_in.(REM_tvec)(this_idx);
                    this_s = (ones(size(A_in.(REM_tvec)(this_idx)))*ii);
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
    
    imagesc(A_in.(REM_tvec), 1:size(act_array,2),  act_array')
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
            
            plot(A_in.(REM_tvec), ((A_in.(sleep_phase)(:,iC)./max(A_in.(sleep_phase)(:,iC)))*.8)+off_set,  'color', c_ord(find(iA == plot_idx),:), 'linewidth', .5)
            
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
    
%     ax(3) = subplot(7,1,6:7);
%     cla
%     hold on
%     leg_val = [];
%     for iA = plot_idx
%         plot(A_in.(REM_tvec), log10(A_in.(REM_proj)(iA,:)), 'color', c_ord(find(iA == plot_idx),:), 'linewidth', 2)
%         
%         leg_val{find(iA == plot_idx)} = ['Assembly #' num2str(iA)];
%     end
%    ylim([0 inf])
% ylabel({'assembly strength'})
% xlabel('time (s)')
% y_val = get(gca, 'YTick'); 
% set(gca, 'YTickLabel', 10.^y_val, 'linewidth', 1); 
% 
% if strcmpi(REM_proj, 'REM_proj')
% 
%             yline(log10(A_in.pREM_stats.R_thresh), '--', 'color', [.7 .7 .7]); 
% 
% else
%         yline(log10(A_in.([REM{1} 'stats']).R_thresh), '--', 'color', [.7 .7 .7]); 
% end
%     legend([leg_val 'Sig thresh.'], 'Orientation', 'horizontal', 'box', 'off')
%     
%     

    %%
    ax(3) = subplot(3,1,3);
cla
hold on
leg_val = [];

for iA = plot_idx
    this_proj = A_in.(REM_proj)(iA,:);
    %         this_proj(this_proj <=0) =NaN;
    plot(A_in.(REM_tvec), log10(this_proj), 'color', c_ord(find(iA == plot_idx),:), 'linewidth', 2)
    leg_val{find(iA == plot_idx)} = [ num2str(iA)];
end
ylim([0 inf])
ylabel({'assembly strength'})
xlabel('time (s)')
y_val = get(gca, 'YTick'); 
set(gca, 'YTickLabel', 10.^y_val, 'linewidth', 1); 
   ylim([0 inf])
    y_val = get(gca, 'YTick');
            y_val = [0 1 2];
%     yTick= [0 1 2 2.3979]
    set(gca, 'ytick', y_val, 'YTickLabel', 10.^y_val, 'linewidth', 1);
% text('R thresh'

 
if strcmpi(REM_proj, 'pREM_proj')

            yline(log10(A_in.pREM_stats.R_thresh), '--', 'color', [.7 .7 .7]); 

else
        yline(log10(A_in.([REM{1} 'stats']).R_thresh), '--', 'color', [.7 .7 .7]); 
end
    legend([leg_val 'Sig thresh.'], 'Orientation', 'horizontal', 'box', 'off')

    linkaxes(ax,'x'); 
    colormap(ax(1), c_map); 
    colormap(ax(3), c_map); 
    
    xlim([A_in.(REM_tvec)(1) A_in.(REM_tvec)(end)])
    
    cfg.ft_size= ft_size;
cfg.resize = 0;
SetFigure(cfg, gcf)
    
    %% save the figure
    if ~isempty(fig_dir)
        pause(.5)
        saveas(gcf, [fig_dir filesep A_in.info.subject '_' A_in.info.session '_' strrep(num2str(A_in.info.bin), '.', 'p') 's_bin_' sleep_phase '_raster.png']);
       pause(.5)
            set(gcf,'Units','Inches');
            pos = get(gcf,'Position');
            set(gcf, 'Position', [pos(1) pos(2) pos(4) pos(4)])
            set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(4), pos(4)])
            ft_size = 14;

    ax(3) = subplot(3,1,2);
            set(gca, 'YTick', c_step(2:end) - diff(c_step)/2, 'YTickLabel',[])
            for ii = 2:length(c_step)
                text(0, c_step(ii) - diff(c_step(ii-1:ii))/2, c_label{ii-1},'HorizontalAlignment', 'right', 'Color', c_map(ii,:), 'FontSize', ft_size )
            end

        print(gcf,  [fig_dir filesep A_in.info.subject '_' A_in.info.session '_' strrep(num2str(A_in.info.bin), '.', 'p') 's_bin_' sleep_phase '_raster.pdf'], '-dpdf','-r0')

    end
