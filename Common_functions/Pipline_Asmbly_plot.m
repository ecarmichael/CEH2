function Pipline_Asmbly_plot(A_out, fig_dir); 


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

for iB = size(A_out,2):-1:1
    MS_Asmbly_plot(A_out{iB}.P_temp, A_out{iB}.P_pos, A_out{iB}.map, A_out{iB}.P_loc, fig_dir, A_out{iB}.info);
close all;
end

%% summarize the reactivations 

for iB = size(A_out,2):-1:1
    MS_Asmbly_ReAct_plot(A_out{iB}, fig_dir);

end

