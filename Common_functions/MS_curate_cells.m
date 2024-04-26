function [keep_idx, ms] = MS_curate_cells(ms)
%% MS_curate_cells: displays the raw traces, binary, deconv (if present) and the spf of each cell across a recording to determine if it is worth keeping.



%%
c_idx = 1:ms.numNeurons;

figure(9999)
set(gcf, 'Units', 'normalized', 'Position', [0.15 0.25 .75 .5])

subplot(2,3,3)
MS_plot_all_SFPs(ms.SFPs_sharp);


c_ord = MS_linspecer(4);

s_idx = nearest_idx3([0 60], ms.time);
e_idx = nearest_idx3([ms.time(end) - 60 ms.time(end)], ms.time);

disp('Keep (space) | reject (delete)')
ii = 1;

keep_idx = zeros(ms.numNeurons, 1);
for kk = 1:ms.numNeurons
    
    
    subplot(2,3,[1 2 ])
    cla
    
    hold on
    plot(ms.time, (ms.RawTraces(:,ii)),'color', [0.6 0.6 0.6], 'linewidth', 1)
    if isfield(ms, 'deconv')
        plot(ms.time, (ms.denoise(:,ii)),'color', [c_ord(1,:) .5], 'linewidth', 1)
        
        this_deconv = ms.deconv(:,ii);
        this_deconv(this_deconv==0) = NaN;
        plot(ms.time, ((this_deconv*5) - 1/4),'color','k', 'linewidth',1)
        
    end
    ylim([-0.5 1.5])
    
    subplot(2,3,4)
    
    cla
    
    hold on
    plot(ms.time(s_idx(1):s_idx(2)), (ms.RawTraces(s_idx(1):s_idx(2),ii)),'color', [0.6 0.6 0.6], 'linewidth', 1)
    if isfield(ms, 'deconv')
        plot(ms.time(s_idx(1):s_idx(2)), (ms.denoise(s_idx(1):s_idx(2),ii)),'color', [c_ord(1,:) .5], 'linewidth', 1)
        
        this_deconv = ms.deconv(s_idx(1):s_idx(2),ii);
        this_deconv(this_deconv==0) = NaN;
        plot(ms.time(s_idx(1):s_idx(2)), ((this_deconv*5) - 1/4),'color','k', 'linewidth',1)
        
    end
    xlim([ms.time(s_idx(1)) ms.time(s_idx(2))])
    ylim([-0.5 1.5])
    
    subplot(2,3,5)
    
    cla
    
    hold on
    plot(ms.time(e_idx(1):e_idx(2)), (ms.RawTraces(e_idx(1):e_idx(2),ii)),'color', [0.6 0.6 0.6], 'linewidth', 1)
    if isfield(ms, 'deconv')
        plot(ms.time(e_idx(1):e_idx(2)), (ms.denoise(e_idx(1):e_idx(2),ii)),'color', [c_ord(1,:) .5], 'linewidth', 1)
        
        this_deconv = ms.deconv(e_idx(1):e_idx(2),ii);
        this_deconv(this_deconv==0) = NaN;
        plot(ms.time(e_idx(1):e_idx(2)), ((this_deconv*5) - 1/4),'color','k', 'linewidth',1)
        
    end
    xlim([ms.time(e_idx(1)) ms.time(e_idx(2))])
    ylim([-0.5 1.5])
    
    
    subplot(2,3,3)
    cla
    imagesc(ms.CorrProj)
    hold on
    [row, col] = find(ms.SFPs_sharp(:,:,ii) == max(max(ms.SFPs_sharp(:,:,ii))));
    t =  scatter(col, row,50,'o','MarkerEdgeColor', 'g');
    %         scatter(col, row,50,'o',  'MarkerEdgeColor',c_ord(ii,:) , 'LineWidth', 1);% c_ord(ii,:)
    title(['Cell: ' num2str(ii)  '/' num2str(ms.numNeurons)]);
    
    
    subplot(2,3,6)
    cla
    imagesc(ms.SFPs_sharp(:,:,ii))
    axis xy
    
    % delete(t)
    
    
    drawnow
    was_a_key = waitforbuttonpress;
    key_hit = get(gcf, 'CurrentKey');
    disp(key_hit)
    if strcmp(key_hit, 'delete') || strcmp(key_hit, 'backspace')
        keep_idx(ii) = 0;
        fprintf('Cell <strong>%0.0f</strong>/%0.0f : %s\n', ii, ms.numNeurons, 'rejected')
        ii = ii+1;
        
        
    elseif strcmp(key_hit, 'return') || strcmp(key_hit, 'space')
        keep_idx(ii) = 1;
        fprintf('Cell <strong>%0.0f</strong>/%0.0f : %s\n', ii, ms.numNeurons, 'accepted')
        ii = ii+1;
        
    elseif strcmp(key_hit, 'leftarrow')
        fprintf('Redo Cell <strong>%0.0f</strong>/%0.0f : %s\n', ii)
                ii = ii-1;

        continue
    else
        
        
    end
    
    
    clear t
end

ms.keep_idx = keep_idx;
end