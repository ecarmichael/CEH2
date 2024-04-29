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


%% find the assemblies with the most spatial tuning and frequent activity

p_rank = []; 

for iB = size(A_out,2):-1:1
    p_cent = []; p_rate = [];

    for ii = length(A_out{iB}.P_pos):-1:1
        p_cent(ii) = A_out{iB}.map{ii}.cent_z;
        p_rate(ii) = length(A_out{iB}.P_loc{ii}.loc_time);
    end
    
    p_rate(p_cent > -1.95) = 0; 
    [~, s_idx] = sort(p_rate, 'descend');
    
    p_rank{iB} = s_idx; 
    
end
%% Plot the raw/raster in time



for iB = size(A_out,2):-1:1
    if length(A_out{iB}.P_pos) >3
        p_idx = p_rank{iB}(1:4);
    elseif (isempty(A_out{iB}.P_pos) || ~isfield(A_out{iB}, 'P_pos'))
        continue
    else
        p_idx = p_rank{iB}(1:length(A_out{iB}.P_pos));
    end
    
    MS_Asmbly_plot_raster(A_out{iB}, fig_dir,p_idx)
    
    MS_Asmbly_plot_raster_figure(A_out{iB}, fig_dir,p_idx); 
    
    
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