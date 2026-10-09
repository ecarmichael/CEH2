function [h, corr_mat] = MS_asmbly_ephys_raster(S, tvec, A_temp, A_proj, idx, tsd_in)
%% :
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
% EC 2025-09-30   initial version 
%
%
%
%% initialize

if nargin < 5
    idx = []; 
end

if isempty(idx)
    idx = 1:size(A_temp,2); 
end

c_ord = MS_linspecer(length(idx));

    %% compute the xcorr between pairs of member and non member cells
    corr_mat = []; 
    cfg.binsize = .01; 
    cfg.max_t = .2; 

    % Snip the data for speed;

    S_snip = restrict(S, tvec(1), tvec(1)+300); %just use 5mins
    for jj = length(S_snip.t):-1:1
        for kk = length(S_snip.t):-1:1
                    fprintf([num2str(jj) '-' num2str(kk)])
            if jj == kk
                corr_mat{jj, kk} = nan;
            else
                tic
                [corr_mat{jj, kk}, t_vec] = ccf(cfg, S_snip.t{jj}, S_snip.t{kk});
                toc
            end
        end
    end

%% plot the assembly templates
figure(909)
clf

m  =  length(idx); 
n =  6;
s_idx = reshape(1:n*m, n, m)'; 
A_cells = []; A_cells_id = []; 
for ii = 1:length(idx)
    subplot(m,n, s_idx(ii,1))
    cla
    hold on
    stem(A_temp(:,idx(ii)), 'color', [.8 .8 .8 .2])

    if isfield()
    a_idx = sum(zscore(A_temp(:,idx(ii))) > 1, 2) > 0;

    stem(find(a_idx), A_temp(find(a_idx),idx(ii)), 'color',c_ord(ii,:), 'MarkerFaceColor', c_ord(ii,:))
    ylim([-.4 .6])
    view(90,90)
    text(0, .6,  ['A: ' num2str(idx(ii))], 'color', c_ord(ii,:), 'VerticalAlignment','top', 'HorizontalAlignment','left')

    A_cells = [A_cells, find(a_idx)']; 
    A_cells_id = [A_cells_id repmat(ii,1,  length(find(a_idx)))];

    subplot(m,n, s_idx(ii,2))
    cla
        hold on

        a = find(a_idx); 
        pairs{1} = 'start'; 
        a_corr = []; 
        for jj = 1:length(a)
            this_corr = [];
            for kk = 1:length(a)
                if jj == kk || sum(contains(pairs,[num2str(kk) '_' num2str(jj)] )) >0
                    continue
                else
                    % disp([num2str(jj) ' ' num2str(kk)])
                    this_corr(end+1,:) =  corr_mat{jj, kk};
                    pairs{kk} = [num2str(jj) '_' num2str(kk)];
                end
            end
            a_corr{jj} = this_corr;
        end

        a_corr_mat = []; 
        for jj = 1:length(a_corr)
            a_corr_mat = [a_corr_mat; a_corr{jj}];
            plot(t_vec, mean(a_corr{jj},1), 'Color',[c_ord(ii,:) .2], 'LineWidth',.1);
        end

        plot(t_vec, mean(a_corr_mat), 'Color',c_ord(ii,:), 'LineWidth',3);

        % for the non members
         nm = find(~a_idx); 
         pairs = []; 
        pairs{1} = 'start'; 

        for jj = 1:length(nm)
            this_corr = [];
            for kk = 1:length(nm)
                if jj == kk || sum(contains(pairs,[num2str(kk) '_' num2str(jj)] )) >0
                    continue
                else
                    disp([num2str(jj) ' ' num2str(kk)])
                    this_corr(end+1,:) =  corr_mat{jj, kk};
                    pairs{kk} = [num2str(jj) '_' num2str(kk)];
                end
            end
            n_corr{jj} = this_corr;
        end

        n_corr_mat = [];
        for jj = 1:length(n_corr)
            n_corr_mat = [n_corr_mat; n_corr{jj}];
            plot(t_vec, mean(n_corr{jj},1), 'Color',[.8 .8 .8 .1], 'LineWidth',.1);
        end
        ylim([0 .2])
end




%% colour code the spikes

cell_ids = 1:length(A_temp); 

rm_idx = ismember(cell_ids,sort(A_cells));

c_ids = [A_cells, cell_ids(~rm_idx)]; 

A_cells_id = [A_cells_id repmat(99,1,  length(find(cell_ids(~rm_idx))))];

s_c_ord = []; 

for ii = 1:length(A_cells_id)
    if A_cells_id(ii) == 99
        s_c_ord(ii,:) = [.7 .7 .7];
    else
        s_c_ord(ii,:) = c_ord(A_cells_id(ii),:);
    end
end

S_sort  = S; 

S_sort.t = S_sort.t(c_ids);
S_sort.label = S_sort.label(c_ids);

%% plot the raster 

ax(1) = subplot(m,n,sort(reshape(s_idx(1:end-1,3:end), 1, numel(s_idx(1:end-1,3:end)))));
cla
cfg_mr.spkColor = s_c_ord; 
cfg_mr.openNewFig = 0; 
h = MultiRaster(cfg_mr, S_sort);


% for iA  = 1:size(A_temp,2)


    % this_S = SelectTS([], S, )


%  h = cell(size(S_sort.t));
%     for iC = 1:nCells
%         nSpikes = length(S.t{iC});
%         xvals = [ S_in.t{iC}' ; S_in.t{iC}' ; nan(size(S_in.t{iC}')) ];
%         yvals = [ iC-cfg.SpikeHeight.*ones(1,nSpikes);...
%             iC+cfg.SpikeHeight.*ones(1,nSpikes) ; nan(1,nSpikes) ];
% 
%         xvals = xvals(:);
%         yvals = yvals(:);
% 
%         if isempty(xvals) && isempty(yvals)
% %             h(iC) = nan;
%             h{iC} = nan;
%         else
% %             h(iC) = plot(xvals,yvals,'Color',cfg.spkColor(iC,:),'LineWidth',cfg.LineWidth);
%             plot(xvals,yvals,'Color',cfg.spkColor(iC,:),'LineWidth',cfg.LineWidth);
%             h{iC} = get(gca);
%         end
%     end




% add some LFP if exists
if exist('tsd_in', 'var')

hold on
data_norm = MS_norm_range(tsd_in.data(1,:), -5, 5); 

    plot(tsd_in.tvec, data_norm-5, 'color', 'k'); 

    ylim([-10 inf])
end


% plot the projections
ax(2) = subplot(m,n,s_idx(end,3:end));
cla
hold on
for ii = 1:length(idx)

    plot(tvec', A_proj(idx(ii),:), 'color', c_ord(ii,:), 'LineWidth',1.5)
end
% set(gca, 'YScale', 'log')
linkaxes(ax, 'x')
xlim([tvec(1) tvec(end)])


%% 