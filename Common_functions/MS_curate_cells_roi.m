function [keep_idx, ms, roi] = MS_curate_cells(ms, ms2, roi)
%% MS_curate_cells: displays the raw traces, binary, deconv (if present) and the spf of each cell across a recording to determine if it is worth keeping.



%% init

if nargin <2
    ms2 =[];
    roi = [];
elseif nargin < 3
    roi = [];
end


%% if an ROI is specified use if to exclude cells outside the ROI. if ROI is just 1 then let the user make one,.
if roi == 1
    figure(101)
    clf
    MS_plot_all_SFPs(ms.SFPs_sharp);
    caxis([0 100])
    poly_roi = MS_drawpoly_wait('', MS_linspecer(1));
    
    roi = [];
    roi.Position = [poly_roi.Position(:,1), poly_roi.Position(:,2)]; % keep the roi positions. If the figure is closed the roi goes too.
end

if ~isempty(roi)
    
    % don't use centroids. instead compute from sharp sfps.
    ms.sharp_centroids = [];
    for ii = size(ms.SFPs_sharp,3):-1:1
        [ms.sharp_centroids(ii,1), ms.sharp_centroids(ii,2)] = find(ms.SFPs_sharp(:,:,ii) == max(max(ms.SFPs_sharp(:,:,ii))));
    end
    
    [c_in, ~] = inpolygon(ms.sharp_centroids(:,2), ms.sharp_centroids(:,1), roi.Position(:,1), roi.Position(:,2));
    hold on
    roi_keep_idx = c_in; % keep cells with centroids in or on the lines;
    scatter(ms.sharp_centroids(roi_keep_idx,2), ms.sharp_centroids(roi_keep_idx,1),25, 'r')
else
    roi_keep_idx = true(1, ms.numNeurons);
end

%% actually check the cells


% c_idx = 1:ms.numNeurons;

figure(9999)
set(gcf, 'Units', 'normalized', 'Position', [0.15 0.25 .75 .5])

subplot(2,3,3)
temp_sfps = ms.SFPs_sharp;

temp_sfps(:,:,~roi_keep_idx) = [];
MS_plot_all_SFPs(temp_sfps);
caxis([0 100])

c_ord = MS_linspecer(4);

s_idx = nearest_idx3([0 60], ms.time);
if isempty(ms2)
    e_idx = nearest_idx3([ms.time(end) - 60 ms.time(end)], ms.time);
else
    e_idx = nearest_idx3([ms2.time(end) - 60 ms2.time(end)], ms2.time);
end

disp('Keep (space) | reject (delete)')
redo = 0;

keep_idx = zeros(ms.numNeurons, 1);
cell_idx = find(roi_keep_idx);

%%
kk = 1; 
while kk <length(cell_idx)
    
    if redo == 0
        kk = kk+1;
    else
        redo =0;
        kk = kk-1;
        if kk <1
            kk = 1;
        end
    end
            ii = cell_idx(kk);

    if isempty(ms2)
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
        
    else
        subplot(2,3,[1 ])
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
        
        subplot(2,3,[2])
        cla
        
        hold on
        plot(ms2.time, (ms2.RawTraces(:,ii)),'color', [0.6 0.6 0.6], 'linewidth', 1)
        if isfield(ms2, 'deconv')
            plot(ms2.time, (ms2.denoise(:,ii)),'color', [c_ord(1,:) .5], 'linewidth', 1)
            
            this_deconv = ms2.deconv(:,ii);
            this_deconv(this_deconv==0) = NaN;
            plot(ms2.time, ((this_deconv*5) - 1/4),'color','k', 'linewidth',1)
            
        end
        ylim([-0.5 1.5])
    end
    
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
    
    
    if isempty(ms2)
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
        
    else
        subplot(2,3,5)
        cla
        
        hold on
        plot(ms2.time(e_idx(1):e_idx(2)), (ms2.RawTraces(e_idx(1):e_idx(2),ii)),'color', [0.6 0.6 0.6], 'linewidth', 1)
        if isfield(ms2, 'deconv')
            plot(ms2.time(e_idx(1):e_idx(2)), (ms2.denoise(e_idx(1):e_idx(2),ii)),'color', [c_ord(1,:) .5], 'linewidth', 1)
            
            this_deconv = ms2.deconv(e_idx(1):e_idx(2),ii);
            this_deconv(this_deconv==0) = NaN;
            plot(ms2.time(e_idx(1):e_idx(2)), ((this_deconv*5) - 1/4),'color','k', 'linewidth',1)
            
        end
        xlim([ms2.time(e_idx(1)) ms2.time(e_idx(2))])
        ylim([-0.5 1.5])
        
        
    end
    
    
    subplot(2,3,3)
    cla
    imagesc(ms.CorrProj)
    hold on
    [row, col] = find(ms.SFPs_sharp(:,:,ii) == max(max(ms.SFPs_sharp(:,:,ii))));
    t =  scatter(col, row,50,'o','MarkerEdgeColor', 'k');
    %         scatter(col, row,50,'o',  'MarkerEdgeColor',c_ord(ii,:) , 'LineWidth', 1);% c_ord(ii,:)
    title(['Cell: ' num2str(ii)  '/' num2str(ms.numNeurons)]);
    
    
    subplot(2,3,6)
    cla
    imagesc(ms.SFPs_sharp(:,:,ii))
    axis xy
    
    % delete(t)
    
    
    drawnow
    
    %         if roi_keep_idx(kk) == 0
    %             keep_idx(ii) = 0;
    %             fprintf('Cell <strong>%0.0f</strong>/%0.0f : %s\n', ii, ms.numNeurons, 'rejected based on ROI')
    %             ii = ii+1;
    %
    %
    %             continue
    %         else
    was_a_key = waitforbuttonpress;
    key_hit = get(gcf, 'CurrentKey');
    disp(key_hit)
    
    if strcmp(key_hit, 'delete') || strcmp(key_hit, 'backspace')
        keep_idx(ii) = 0;
        fprintf('Cell <strong>%0.0f</strong>/%0.0f : %s\n', ii, ms.numNeurons, 'rejected')
        %                 ii = ii+1;
        
        
    elseif strcmp(key_hit, 'return') || strcmp(key_hit, 'space')
        keep_idx(ii) = 1;
        fprintf('Cell <strong>%0.0f</strong>/%0.0f : %s\n', ii, ms.numNeurons, 'accepted')
        %                 ii = ii+1;
        
    elseif strcmp(key_hit, 'leftarrow')
        fprintf('Redo Cell <strong>%0.0f</strong>/%0.0f : %s\n', ii)
        %ii = cell_idx(kk-1);
        redo = 1;
        continue
    end
    
    
    clear t
end

ms.keep_idx = keep_idx;
end