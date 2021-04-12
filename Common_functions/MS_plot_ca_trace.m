function MS_plot_ca_trace(data_in, tvec, threshold, linesize)
%% MS_plot_ca_trace:
%
%
%
%    Inputs: 
%    - data_in: [nChan x nSample] matrix with binary events.  
%
%    - tvec: [1 x nSample] time vector [optional]
%
%
%
%    Outputs: 
%    - h [handle] 
%
%
% EC 2021-03-27   initial version 
%
%% initialize

if nargin < 2 || isempty(tvec)
    tvec = 1:size(data_in,2);
    fprintf('No tvec give, using data length (%i samples) as vector\n', size(data_in, 2))
end

if nargin <3
    threshold = 0;
end

if nargin < 4
    linesize = 2; 
end

if size(tvec,2) ~= size(data_in, 2)
    error('data_in and tvec differ in size')
end
%%
nChan = size(data_in, 1); 
% set colours
c_ord = linspecer(nChan+1); % add one more for background
j_ord = jet(100); % set the set 



% set a threshold for plotting.  for clearity. 
data_in(data_in <= threshold) = NaN; 

% figure
hold on
% set(gca, 'YDir','reverse'); % flip y axis order
xlim([tvec(1) tvec(end)]);
ylim([1, nChan+1])

for iC = nChan:-1:1
%     keepIndex = ~isnan(data_in(iC,:)+iC);

%     ax = patch([tvec , fliplr(tvec)], [data_in(iC,:)+iC,(iC-1)*ones(size(data_in(iC,:)))], c_ord(iC,:)); 
%     ax.EdgeColor = [0 0 0]; 
%     ax = area(data_in(iC,:)+iC, iC, 'ShowBaseLine', 'off'); 
%     ax.EdgeColor = c_ord(iC,:);
%     ax.FaceColor = c_ord(iC,:);

    plot(tvec, (data_in(iC,:)./max(data_in(iC,:)))+iC+2, 'color', c_ord(iC,:), 'LineWidth', linesize)
    tick_val(iC) = nanmedian((data_in(iC,:)/max(data_in(iC,:)))+iC+1);
    tick_label{iC} = num2str(iC);
end

    
%     
%     
% %% set up defaults.
% cfg_def =[];
% cfg_def.Ca_type = 'RawTraces'; % can be either 'RawTraces' or FiltTraces' (maybe others?)
% cfg_def.Ca_chan = 1:floor(n_cells/20):n_cells; % get a subset of cells.
% cfg_def.plot_type = '3d'; % '2d' or '3d'
% cfg_def.x_zoom = []; % where to zoom in on the x_axis for each plot.
% cfg_def.view = [0 45]; 
% cfg_def.colors = []; % can be used to pass a colour set. 
% cfg_def.label = 'file_names';
% cfg_def.saveas = []; % if empty don't save the images.  Can be '.fig' (matlab fig) or other known saveas format.  I like '.png'.
% 
% cfg = ProcessConfig(cfg_def, cfg_in);
% 
% if n_cells < 30
%     cfg.Ca_chan = 1:n_cells;
% end
% 
% %% make the plots
% if isempty(cfg.colors)
%     c_basic = linspecer(10); % basic color range
%     c_basic(5:6,:) = [];
%     c_ord = repmat(c_basic, (ceil(size(cfg.Ca_chan,2)/length(c_basic))),1); % repeat the colors range for better visibility;
%     % c_ord = linspecer(length(cfg.Ca_chan)); % nice colours.
% else
%     c_ord = cfg.colors;
% end

        

% 
% xlim([time_in2(1)*0.001 time_in2(end)*0.001]);
% y_lim = ylim;
% ylim([0 y_lim(2)]); 
% 
% %         ax = [];
% if ~isempty(cfg.x_zoom)
%     
%     xlim(cfg.x_zoom);
% end
% 
% switch cfg.plot_type
%     case '2d'
%         [~, sort_idx] = sort(tick_val); % get past cases where the signal is lower than the offset previous. 
%         set(gca, 'ytick', tick_val(sort_idx), 'yticklabel', tick_label(sort_idx));
%         ylabel('cell id')
end % end function