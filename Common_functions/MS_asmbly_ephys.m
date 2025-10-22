function out = MS_asmbly_ephys(data_in, move_idx, opts)
%% MS_asmbly_ephys: provides a wrapper for running assembly and reactivation analyses using Ts spike data data.

if nargin <1
    error('needs spike and tvec')
elseif nargin < 2
    opts = [];
elseif nargin < 3
    opts = []; 
end



if isnumeric(data_in)
    fprintf('<strong>%s</strong>: data in is numeric, treating as if it is already binned\n', mfilename)
    data_h = data_in;

elseif isfield(data_in, 'type')
    if strcmpi(data_in.type, 'tsd')
       fprintf('<strong>%s</strong>: data in is TSD, treating as if it is already binned\n', mfilename)
       data_h = data_in.data; 
    end
end


if isempty(opts)
   opts.threshold.method = 'MarcenkoPastur';
   % opts.threshold.method = 'circularshift';
    opts.Patterns.method = 'ICA';
    opts.Patterns.number_of_iterations = 500;
    opts.threshold.number_of_permutations = 500; 
    opts.threshold.permutations_percentile = 95;
end

rng(123, 'twister')
%% convert the TS data into rates

%% try the assembly code....


% get the assemblies

A_temp = assembly_patterns(data_h', opts);

if size(A_temp,1) == size(data_h',2)
    disp('Templates match # cells')
else 
    error('Input data has more cells than timepoints');
end


% get the projections
rng(123, 'twister')
A_proj = assembly_activity(A_temp,data_h);


%% shuffle distribution for assemblies
rng(123,'twister')
nShuff = 100;
wake_shuff_mat = [];

Ass_shuff = NaN(1,nShuff);
for iS = nShuff:-1:1
    tic
    shuff_data = NaN(size(data_h));
    for ic = 1:size(data_h,2)
        shuff_data(:,ic) = circshift(data_h(:,ic), floor(MS_randn_range(1,1,1,size(data_h,1))));
    end
    
    this_ass = assembly_patterns(shuff_data');
    if ~isempty(this_ass)
        S_prog = assembly_activity(this_ass,shuff_data');
        
        wake_shuff_mat(iS,:) =  S_prog(1,:);
        keep_idx(iS) = 1;
    else
        wake_shuff_mat(iS,:) = NaN(1,length(shuff_data));
        keep_idx(iS) = 0;
    end
    %     for ii = size(this_ass,2):-1:1
    
    if sum(max(this_ass) > 0.2) >0
        Ass_shuff(iS) = sum(max(this_ass) > 0.2);
    else
        Ass_shuff(iS) = 0;
    end
    %     end
    % fprintf('Shuff # %.0f found %.0f assemblies and took %2.2f seconds\n', iS, size(this_ass,2), toc)
end

shuff_stats.shuff_n = Ass_shuff; 
shuff_stats.mean = mean(Ass_shuff); 
shuff_stats.sd = std(Ass_shuff); 
shuff_stats.p95 = prctile(Ass_shuff, 95, 'all');
shuff_stats.p99 = prctile(Ass_shuff, 99, 'all');



w_thresh = prctile(wake_shuff_mat(wake_shuff_mat >0), 99, 'all');

