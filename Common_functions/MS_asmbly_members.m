function [asmbly_out] = MS_asmbly_members(asmbly_in, method)
%% get the cells with significant weights per assembly

if nargin < 2
    method = 'SD'; 
end
asmbly_out = asmbly_in;
asmbly_out.P_cells = P_cells;
asmbly_out.P_cells_id = P_cells_id;

fprintf('Finding members using: <strong>%s</strong>\n', method)

idx = 1:size(asmbly_in.P_temp,2); 

% set the colours
c_ord = MS_linspecer(length(idx)); 

m  =  ceil(length(idx)/4);
n  =  ceil(length(idx)/4);
s_idx = 1:(n*m); 
P_cells = []; P_cells_id = []; 

for ii = 1:length(idx)
    subplot(m,n, s_idx(ii))
    cla
    hold on
    stem(asmbly_in.P_temp(:,idx(ii)), 'color', [.8 .8 .8 .2])

    if strcmpi(method, 'otsu')
        % [o_data o_idx] = sort(abs(asmbly_in.P_temp(:,idx(ii)))); 
        [o_data] = (abs(asmbly_in.P_temp(:,idx(ii)))); 
        % o_data(o_data<0) = NaN; 
        o_thresh = otsuthresh(o_data); 
        % convert back to data scale
        o_thresh = o_thresh*(max(o_data)-min(o_data)) + min(o_data); 
        a_idx = sum(asmbly_in.P_temp(:,idx(ii)) > o_thresh, 2) > 0;

    else
    a_idx = sum(zscore(asmbly_in.P_temp(:,idx(ii))) > 1, 2) > 0;

    end

    stem(find(a_idx), asmbly_in.P_temp(find(a_idx),idx(ii)), 'color',c_ord(ii,:), 'MarkerFaceColor', c_ord(ii,:))
    ylim([-.4 .6])
    view(90,90)
    text(0, .6,  ['A: ' num2str(idx(ii))], 'color', c_ord(ii,:), 'VerticalAlignment','top', 'HorizontalAlignment','left')

    P_cells = [P_cells, find(a_idx)']; 
    P_cells_id = [P_cells_id repmat(ii,1,  length(find(a_idx)))];

    % subplot(m,n, s_idx(ii,2))
    % cla
    %     hold on
    % 
    %     a = find(a_idx); 
    %     pairs{1} = 'start'; 
    %     a_corr = []; 
    %     for jj = 1:length(a)
    %         this_corr = [];
    %         for kk = 1:length(a)
    %             if jj == kk || sum(contains(pairs,[num2str(kk) '_' num2str(jj)] )) >0
    %                 continue
    %             else
    %                 % disp([num2str(jj) ' ' num2str(kk)])
    %                 this_corr(end+1,:) =  corr_mat{jj, kk};
    %                 pairs{kk} = [num2str(jj) '_' num2str(kk)];
    %             end
    %         end
    %         a_corr{jj} = this_corr;
    %     end
    % 
    %     a_corr_mat = []; 
    %     for jj = 1:length(a_corr)
    %         a_corr_mat = [a_corr_mat; a_corr{jj}];
    %         plot(t_vec, mean(a_corr{jj},1), 'Color',[c_ord(ii,:) .2], 'LineWidth',.1);
    %     end
    % 
    %     plot(t_vec, mean(a_corr_mat), 'Color',c_ord(ii,:), 'LineWidth',3);
    % 
    %     % for the non members
    %      nm = find(~a_idx); 
    %      pairs = []; 
    %     pairs{1} = 'start'; 
    % 
    %     for jj = 1:length(nm)
    %         this_corr = [];
    %         for kk = 1:length(nm)
    %             if jj == kk || sum(contains(pairs,[num2str(kk) '_' num2str(jj)] )) >0
    %                 continue
    %             else
    %                 disp([num2str(jj) ' ' num2str(kk)])
    %                 this_corr(end+1,:) =  corr_mat{jj, kk};
    %                 pairs{kk} = [num2str(jj) '_' num2str(kk)];
    %             end
    %         end
    %         n_corr{jj} = this_corr;
    %     end
    % 
    %     n_corr_mat = [];
    %     for jj = 1:length(n_corr)
    %         n_corr_mat = [n_corr_mat; n_corr{jj}];
    %         plot(t_vec, mean(n_corr{jj},1), 'Color',[.8 .8 .8 .1], 'LineWidth',.1);
    %     end
    %     ylim([0 .2])
end