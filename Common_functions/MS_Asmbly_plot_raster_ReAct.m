function MS_Asmbly_plot_raster_ReAct(A_in, fig_dir,sleep_phase,plot_idx, type)

if nargin < 3
    fig_dir = [];
    plot_idx = 1:length(A_in.P_pos);
    type = 'act_mat';
elseif nargin < 4
    plot_idx = 1:length(A_in.P_pos);
    type = 'act_mat';
elseif nargin < 5
    type = 'act_mat';
end


REM = strsplit(sleep_phase, 'data'); 
    REM_tvec = [ REM{1} 'tvec'];
    REM_proj = [REM{1} 'proj']; 
%%

c_ord = MS_linspecer(length(plot_idx)+ceil(length(plot_idx)/5));

% make a colour coded activity array
act_array = zeros(size(A_in.(sleep_phase)));

cnt = 0; c_list = []; c_map = [0 0 0]; c_step = 0;
for iA = plot_idx
    c_idx = A_in.P_pos{iA};
    A_range = 1:length(c_idx);
    A_range = A_range+cnt;
    cnt = cnt+length(A_range)+1;
    
    for ii = 1:length(c_idx)
        act_array(logical(A_in.(sleep_phase)(:,c_idx(ii))),A_range(ii)) = find(iA ==plot_idx);
    end
    
    
    c_list = [c_list; c_idx];
    c_map = [c_map; c_ord(find(iA ==plot_idx),:)];
    c_step = [c_step cnt];
    c_label{find(iA ==plot_idx)} = ['A# ' num2str(iA)];
end


%%


figure(900)
clf
maximize

if isstr(type) == 1
    ax(1) = subplot(7,1,1:5);
    cla
    hold on
    off_set = 0; these_idx = [];
    
    if strcmpi(type, 'act_mat')
        imagesc(A_in.(REM_tvec), 1:size(act_array,2),  act_array')
        ylim([0 length(c_list)])
        set(gca, 'YTick', c_step(2:end) - diff(c_step)/2, 'YTickLabel',c_label)
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
    end
    
else
    ax(1) = subplot(7,1,1:3);
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
                
                    plot(A_in.(REM_tvec), ((A_in.(sleep_phase)(:,iC)./max(A_in.wake_data(:,iC)))*.8)+off_set,  'color', c_ord(find(iA == plot_idx),:), 'linewidth', .5)
            
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
    
    ax(3) = subplot(7,1,6:7);
    cla
    hold on
    leg_val = [];
    for iA = plot_idx
        plot(A_in.(REM_tvec), A_in.(REM_proj)(iA,:), 'color', c_ord(find(iA == plot_idx),:), 'linewidth', 2)
        
        leg_val{find(iA == plot_idx)} = ['Assembly #' num2str(iA)];
    end
    ylim([0 inf])
    ylabel({'assembly strength'})
    xlabel('time (s)')
    
    legend(leg_val, 'Orientation', 'horizontal', 'box', 'off')
    
    yline(A_in.([REM{1} 'stats']).R_thresh); 
    
    linkprop(ax,{'XLim'}); 
    colormap(ax(1), c_map); 
    colormap(ax(3),c_map); 
    
    xlim([A_in.(REM_tvec)(1) A_in.(REM_tvec)(end)])
    %% save the figure
    if ~isempty(fig_dir)
        saveas(gcf, [fig_dir filesep A_in.info.subject '_' A_in.info.session '_' strrep(num2str(A_in.info.bin), '.', 'p') 's_bin_' sleep_phase '_raster.png']);
    end
