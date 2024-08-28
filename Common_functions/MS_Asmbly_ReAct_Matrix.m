function ReAct_out = MS_Asmbly_ReAct_Matrix(this_data, plot_flag)
%% MS_Asmbly_ReAct_Matrix: count the number of 'siginificant' replays per training - testing pairs (PRE, Wake, Post).
%
%
%
%    Inputs:
%    - Asmbly struct with both pREM and postREM appended.
%
%
%
%    Outputs:
%    -
%
%
%
%
% EC 2024-08-26   initial version
%
%
%% initialize

if nargin < 2
    plot_flag = 0;
end
%% split the data into target and reference sets
if isempty(this_data.pREM_temp)
    
    
ReAct_out = [];
ReAct_out.nReAct = NaN(3,3);
ReAct_out.nA_ReAct =  NaN(3,3);
ReAct_out.nReAct_Rate =  NaN(3,3);
ReAct_out.A_ReAct_Rate =  NaN(3,3);
ReAct_out.info.labels = {'Post_Post','Post_Pre','Post_Wake';...
    'Pre_Post','Pre_Pre','Pre_Wake';...
    'Wake_Post','Wake_Pre','Wake_Wake'};
ReAct_out.info.descriptions{1} = 'nReAct: number of reactivations [sum(this_nReact(this_nReact >0))]';
ReAct_out.info.descriptions{2} = 'nA_ReAct: number of assemblies reactivated [length(this_nReact(this_nReact >0))]';
ReAct_out.info.descriptions{3} = 'nReAct_Rate: Rate(nA ./ length of recording) [sum(this_nReact(this_nReact >0))./ (length(data.(states{ii}).(states{jj}).proj)*this_data.bins)]';
ReAct_out.info.descriptions{4} = 'A_ReAct_Rate: ReAct Rate per assemblys [mean(this_nReact(this_nReact >0)./ (length(data.(states{ii}).(states{jj}).proj)*this_data.bins))]';
    
else
% Pre.pos = this_data.pREM_A_pos;
% Pre.temp = this_data.pREM_temp;
thresh.Pre = this_data.pREM_stats.R_thresh;
thresh.Wake.Wake = 10;
thresh.Pre = this_data.pREM_stats.R_thresh;


data.Pre.Pre.proj = this_data.pREM_proj;
data.Pre.Wake.proj = this_data.pREM_Wake_proj;
% get the post proj
data.Pre.Post.proj = assembly_activity(this_data.pREM_temp, this_data.REM_Post_data');


% wake
data.Wake.Pre.proj = this_data.REM_Pre_proj;
data.Wake.Wake.proj = this_data.P_proj;
data.Wake.Post.proj =  this_data.REM_Post_proj;

% post
data.Post.Pre.proj = assembly_activity(this_data.postREM_temp, this_data.REM_Pre_data');
data.Post.Wake.proj = this_data.pREM_Wake_proj;
% get the post proj
data.Post.Post.proj = this_data.postREM_proj;

%% loop over states and make a matrix

states = {'Post','Pre','Wake'};
ReAct_mat = []; ReAct_n = []; ReAct_nA_ep = []; ReAct_rate = [];
labels = {};


for ii = 1:length(states)
    for jj = 1:length(states)
        
        this_nReact = sum(data.(states{ii}).(states{jj}).proj > 10,2);
        
        ReAct_mat(jj, ii) = sum(this_nReact(this_nReact >0));
        
        ReAct_n(jj, ii) = length(this_nReact(this_nReact >0));

        ReAct_pA(jj, ii) = sum(this_nReact >0)./length(this_nReact);

        ReAct_nA_ep(jj, ii) = sum(this_nReact(this_nReact >0))./ (length(data.(states{ii}).(states{jj}).proj)*this_data.bins);
        
        ReAct_rate(jj, ii) = mean(this_nReact(this_nReact >0)./ (length(data.(states{ii}).(states{jj}).proj)*this_data.bins));
        
        labels{jj, ii} = ['tar-' states{jj} ' : ref-' states{ii}];
    end
end


if plot_flag
    figure(767)
    subplot(1,5,1)
    imagesc(ReAct_mat);
    set(gca, 'xtick', 1:3, 'xticklabels', states, 'ytick', 1:3, 'yticklabels', states)
    title('number of reactivations'); colorbar;
    ylabel('Target')
    xlabel('Reference')
    
    subplot(1,5,2)
    imagesc(ReAct_n);
    set(gca, 'xtick', 1:3, 'xticklabels', states, 'ytick', 1:3, 'yticklabels', states)
    title('number of assemblies reactivated'); colorbar;
    
        subplot(1,5,3)
    imagesc(ReAct_pA);
    set(gca, 'xtick', 1:3, 'xticklabels', states, 'ytick', 1:3, 'yticklabels', states)
    title('number of assemblies reactivated'); colorbar;
    
    subplot(1,5,4)
    imagesc(ReAct_nA_ep);
    set(gca, 'xtick', 1:3, 'xticklabels', states, 'ytick', 1:3, 'yticklabels', states)
    title('Rate(nA ./ length of recording)'); colorbar;
    
    subplot(1,5,5)
    imagesc(ReAct_rate);
    set(gca, 'xtick', 1:3, 'xticklabels', states, 'ytick', 1:3, 'yticklabels', states)
    title('ReAct Rate per assembly'); colorbar;
end

ReAct_out = [];
ReAct_out.nReAct = ReAct_mat;
ReAct_out.nA_ReAct = ReAct_n;
ReAct_out.nReAct_Rate = ReAct_nA_ep;
ReAct_out.A_ReAct_Rate = ReAct_rate;
ReAct_out.info.labels = labels;
ReAct_out.info.descriptions{1} = 'nReAct: number of reactivations [sum(this_nReact(this_nReact >0))]';
ReAct_out.info.descriptions{2} = 'nA_ReAct: number of assemblies reactivated [length(this_nReact(this_nReact >0))]';
ReAct_out.info.descriptions{3} = 'nReAct_Rate: Rate(nA ./ length of recording) [sum(this_nReact(this_nReact >0))./ (length(data.(states{ii}).(states{jj}).proj)*this_data.bins)]';
ReAct_out.info.descriptions{4} = 'A_ReAct_Rate: ReAct Rate per assemblys [mean(this_nReact(this_nReact >0)./ (length(data.(states{ii}).(states{jj}).proj)*this_data.bins))]';
end
