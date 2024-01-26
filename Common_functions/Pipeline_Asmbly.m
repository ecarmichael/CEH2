function out = Pipeline_Asmbly(fname,bin_s, move_thresh)
%% Pipeline_Asmbly: provides a wrapper for running assembly and reactivation analyses using calcium data.


if nargin <2
    bin_s = .5;
    move_thresh = 5;
elseif nargin <3
    move_thresh = 5;
end


rng(123, 'twister')


%% load the data

% load('ms_trk.mat')
% load('behav_DLC.mat')
this_sess = fname;

load(this_sess);

%% preprocess beahviour
behav = MS_align_data(behav,ms);

move_idx = behav.speed > move_thresh;

%% remove questionable cells

ms_trk = ms;
keep_idx = zeros(1,size(ms_trk.RawTraces,2));
keep_idx(1:floor(size(ms_trk.RawTraces,2)*.66)) = 1;

remove_cell_id = find(~keep_idx);

cfg_rem = [];
cfg_rem.remove_idx = find(~keep_idx);
cfg_rem.data_type = 'RawTraces';
ms_trk_cut = MS_Remove_trace(cfg_rem, ms_trk);

if ~isfield(ms_trk, 'deconv') % get the deconvolved trace if not already present.
    ms_trk_cut = MS_append_deconv(ms_trk_cut, 1);
end
% remove inactive cells

keep_idx = sum(ms_trk_cut.Binary, 1) >0;

remove_cell_id_decon = find(~keep_idx);

cfg_rem = [];
cfg_rem.remove_idx = find(~keep_idx);
cfg_rem.data_type = 'deconv';
ms_trk_cut = MS_Remove_trace(cfg_rem, ms_trk_cut);


% trim the REM data as well
all_binary_post_REM(:, remove_cell_id) = [];
all_binary_post_REM(:, remove_cell_id_decon) = [];

% all_deconv_post_REM

all_binary_pre_REM(:, remove_cell_id) = [];
all_binary_pre_REM(:, remove_cell_id_decon) = [];

%% load the place information

% dir_parts = strsplit(this_sess, filesep);
% parts = strsplit(dir_parts{end}, '_');
% info.task = parts{2};
% if contains(info.task, 'HATD6')
%     info.task = 'HATDSwitch';
% end
% info.subject = parts{1};


load([info.subject '_' info.session '_PCs.mat'])

place = [];

place.centroids = PCs_properties.peak_loc;

% is it a place cell?
place.is = PCs_properties.isPC;

place.map = PCs_properties.tuning_curve_data';

place.MI = PCs_properties.MI;

place.peak_rate = PCs_properties.peak_rate;



place.centroids(remove_cell_id)= [];
place.centroids(remove_cell_id_decon)= [];

place.is(remove_cell_id)= [];
place.is(remove_cell_id_decon)= [];

place.map(remove_cell_id,:)= [];
place.map(remove_cell_id_decon,:)= [];

place.MI(remove_cell_id)= [];
place.MI(remove_cell_id_decon)= [];

place.peak_rate(remove_cell_id)= [];
place.peak_rate(remove_cell_id_decon)= [];

bin = 3; %80/size(place.map,2);
p_bins = 0:bin:100;
place.p_bins = p_bins(1:end)+bin/2;
% see if there are any anxiety cells

% [~,p_sort] = sort(place.centroids);

%% get the initial assemblies

%     opts.threshold.method = 'binshuffling';
%     opts.Patterns.method = 'ICA';
%     opts.Patterns.number_of_iterations = 500;
%     opts.threshold.permutations_percentile= 95;
%     opts.threshold.number_of_permutations= 500;

A_temp = []; A_prog = []; wake_data = []; wake_tvec = [];
for ii = length(bin_s):-1:1
    
    [A_temp{ii}, A_proj{ii}, wake_data{ii}, wake_tvec{ii}] = MS_PCA_ICA_only(ms_trk_cut, move_idx, bin_s(ii),[]);
end


for ii = length(bin_s):-1:1
    fprintf('PCA-ICA detected %.0f assemblies using a %.2fs binsize\n', size(A_temp{ii},2), bin_s(ii))
end


%% remove assemblies without positive weights
for iB = length(A_temp):-1:1
    
    [P_temp{iB}, P_proj{iB}, P_pos{iB}] = MS_Asmbly_select(A_temp{iB}, A_proj{iB}, 2);
    
    fprintf('[%.0f/%.0f = %.0f%%] Assemblies had cells with positive weights (%.2fs binsize)\n',size(P_temp{iB},2),size(A_temp{iB},2),  (size(P_temp{iB},2)/size(A_temp{iB},2))*100, bin_s(iB))
end


%% get the spacial tuning of the assemblies
min_N_place = 3;

Place_temp = []; Place_proj = []; Place_map = [];
for iB = length(P_temp):-1:1
    
    [map_out{iB}, place_idx{iB}] = MS_Asmbly_map(P_pos{iB}, place, min_N_place);
    
    Place_map{iB} = map_out{iB};
    Place_map{iB}(~place_idx{iB}) = [];
    
    Place_temp{iB} = P_temp{iB}(:,place_idx{iB});
    Place_proj{iB} = P_proj{iB}(:,place_idx{iB});
    
    fprintf('[%.0f/%.0f = %.0f%%] Assemblies contained at least %0.0f place cells (%.2fs binsize)\n',size(Place_temp{iB},2),size(A_temp{iB},2),  (size(Place_temp{iB},2)/size(A_temp{iB},2))*100, min_N_place, bin_s(iB))
    
end

%% get the activation locations on the track;
win_s = 2;
thresh = 10;
for iB = length(P_temp):-1:1
    
    [P_loc{iB}] = MS_Asmbly_act_loc(P_proj{iB}, wake_tvec{iB}, behav, win_s, thresh, 2/bin_s(iB));
    
end


%% Load the REM data and compare it to the wake
cfg_ReAct.nShuff = 500;
cfg_ReAct.thresh = 99;

for iB = length(bin_s):-1:1
    % pre REM
    [REM_pre_proj{iB},REM_pre_stats{iB}, REM_pre_data{iB}, REM_pre_tvec{iB}, REM_Pre_shuff{iB}] = MS_Asmbly_ReAct(cfg_ReAct, all_binary_pre_REM, P_temp{iB} ,ms_trk_cut,  bin_s(iB));
    
    % post REM
    [REM_post_proj{iB}, REM_post_stats{iB}, REM_post_data{iB}, REM_post_tvec{iB},REM_Post_shuff{iB}] = MS_Asmbly_ReAct(cfg_ReAct,all_binary_post_REM, P_temp{iB} ,ms_trk_cut,  bin_s(iB));
    
end


%% get the cross correlation between assemblies
xc_bin = 0.1; 
t_max = 5; 

for iB = length(bin_s):-1:1
    
    [wake_xcor] = MS_Asmbly_xcor(P_proj{iB},wake_tvec, 5, xc_bin, t_max)
    
    
end



%% reactivation strength

for iB = length(bin_s):-1:1
    
    for ii = size(P_proj{iB},1):-1:1
        
        ReAct{iB}(ii) = mean(REM_post_proj{iB}(ii,:)) - mean(REM_pre_proj{iB}(ii,:));
    end
end

%% collect the outputs
for iB = length(bin_s):-1:1
    info.bin = bin_s(iB);
    
    out{iB}.info = info;
    out{iB}.bins = bin_s(iB);
    
    % all templates
    out{iB}.A_temp = A_temp{iB};
    out{iB}.A_proj = A_proj{iB};
    
    % positive templates
    out{iB}.P_temp = P_temp{iB};
    out{iB}.P_proj = P_proj{iB};
    out{iB}.P_pos = P_pos{iB};
    out{iB}.P_loc = P_loc{iB};
    out{iB}.map = map_out{iB};
    
    % place cell assemblies only.
    out{iB}.Place_temp = Place_temp{iB};
    out{iB}.Place_proj = Place_proj{iB};
    out{iB}.Place_map = Place_map{iB};
    
    % reactivations
    out{iB}.ReAct = ReAct{iB};
    
    out{iB}.REM_Pre_proj = REM_pre_proj{iB};
    out{iB}.REM_Pre_stats = REM_pre_stats{iB};
    out{iB}.REM_Pre_data = REM_pre_data{iB};
    out{iB}.REM_Pre_tvec = REM_pre_tvec{iB};
    out{iB}.REM_Pre_shuff = REM_Pre_shuff{iB};
    
    
    out{iB}.REM_Post_proj = REM_post_proj{iB};
    out{iB}.REM_Post_stats = REM_post_stats{iB};
    out{iB}.REM_Post_data = REM_post_data{iB};
    out{iB}.REM_Post_tvec = REM_post_tvec{iB};
    out{iB}.REM_Post_shuff = REM_Post_shuff{iB};
    
end



%% check the correlation and pool



% A_corr = [];
% int_tvec = ms.time(1)/1000:min(bin_s):ms.time(end)/1000;
% int_tvec = int_tvec(1:end-1);
%
% for ii = length(A_prog):-1:1
%
%     A_1_int = []; A_2_int = [];
%
%     for iA = size(A_prog{ii},1):-1:1
%
%         A_1_int = interp1(wake_tvec{ii}, A_prog{ii}(iA,:), int_tvec);
%
%
%
%         for jj = length(A_prog):-1:1
%
%             for iA2 = size(A_prog{jj},1):-1:1
%
%                 A_2_int = interp1(wake_tvec{jj}, A_prog{jj}(iA2,:),int_tvec);
%
%
%                 [this_corr this_p] = corrcoef(A_1_int(~isnan(A_1_int)), A_2_int(~isnan(A_2_int)));
%                 A_corr(iA,iA2) = this_corr(1,2);
%                 P_corr(iA, iA2) = this_p(1,2);
%
%             end
%
%         end
%
%     end
%
% end
%
% sum(P_corr(logical(triu(ones(size(P_corr)),1))) < 0.05)/length(P_corr(logical(triu(ones(size(P_corr)),1))))
