function Pipline_Asmbly_plot_SWS(A_out, fig_dir);


%% initialize

if nargin <2
    save_fig = 0;
elseif nargin == 2
    save_fig = 1;
    if ~exist(fig_dir)
        mkdir(fig_dir)
    end
end



%% plot the summary for all positive cells

%% summarize the reactivations

for iB = size(A_out,2):-1:1
    MS_Asmbly_ReAct_plot(A_out{iB}, fig_dir, 'SWS');
    
end


%% Plot the raw/raster in time

for iB = size(A_out,2):-1:1
    if length(A_out{iB}.P_pos) >3
        p_idx = 1:4;
    elseif (isempty(A_out{iB}.P_pos) || ~isfield(A_out{iB}, 'P_pos'))
        continue
    else
        p_idx = 1:length(A_out{iB}.P_pos);
    end
    
    MS_Asmbly_plot_raster(A_out{iB}, fig_dir,p_idx)
    
    
end


%% plot the raw/raster for the sleep periods
for iB = size(A_out,2):-1:1
    
    if length(A_out{iB}.P_pos) >3
        p_idx = 1:4;
    elseif (isempty(A_out{iB}.P_pos) || ~isfield(A_out{iB}, 'P_pos'))
        continue
    else
        p_idx = 1:length(A_out{iB});
    end
    
    MS_Asmbly_plot_raster_ReAct(A_out{iB},fig_dir,'SWS_Pre_data',p_idx)
    
    
    MS_Asmbly_plot_raster_ReAct(A_out{iB},fig_dir,'SWS_Post_data',p_idx)
    
end