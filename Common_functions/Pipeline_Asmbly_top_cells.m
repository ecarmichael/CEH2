function out = Pipeline_Asmbly_top_cells(fname,bin_s, move_thresh, method)
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

%% load the selected neurons h5
h5_dir = dir('selected*.h5'); 
for ii = length(h5_dir):-1:1
    if  contains(h5_dir(ii).name, this_sess(1:6)) && contains(h5_dir(ii).name,this_sess(8:12))
        keep_idx(ii) = true; 
    else
        keep_idx(ii) = false;
    end

end
h5_idx = find(keep_idx); 
fprintf('session: <strong>%s</strong> found h5 <strong>%s</strong>\n', this_sess, h5_dir(h5_idx).name); 

S_neurons  = h5read(h5_dir(h5_idx).name, '/place_cells'); 
S_neurons = S_neurons+1; 

%% remove questionable cells

ms_trk = ms;

rm_idx = 1:ms_trk.numNeurons; 
remove_cell_id = find(~ismember(rm_idx,double(S_neurons)));


cfg_rem = [];
cfg_rem.remove_idx = remove_cell_id;
cfg_rem.data_type = 'RawTraces';
ms_trk_cut = MS_Remove_trace(cfg_rem, ms_trk);
ms_trk_cut.cell_id = find(ismember(rm_idx,double(S_neurons))); 


% super hack way to sort based on h5 criteria b/c I can't remember how to
% sort this stuff. 
s_idx = []; 
for ii = length(S_neurons):-1:1
   s_idx(ii) = find(S_neurons(ii) == ms_trk_cut.cell_id);
end



if  strcmpi(method, 'grosmark')
    if ~isfield(ms_trk, 'deconv') % get the deconvolved trace if not already present.
        ms_trk_cut = MS_append_deconv(ms_trk_cut, 1);
    end
end

% trim the REM data as well
all_binary_post_REM(:, remove_cell_id) = [];
all_binary_post_REM = all_binary_post_REM(:, s_idx); 


all_binary_pre_REM(:, remove_cell_id) = [];
all_binary_pre_REM = all_binary_pre_REM(:, s_idx); 

% if using 'grosmark 2021' method, create the Csp 
if strcmpi(method, 'grosmark')
    % trim the REM data as well
    all_deconv_pre_REM(:, remove_cell_id) = [];
    all_deconv_pre_REM = all_deconv_pre_REM(:, s_idx);
    
    all_denoise_pre_REM(:, remove_cell_id) = [];
    all_denoise_pre_REM = all_denoise_pre_REM(:, s_idx);
    
    % all_deconv_post_REM
    
    all_deconv_post_REM(:, remove_cell_id) = [];
    all_deconv_post_REM = all_deconv_post_REM(:, s_idx);
    
    all_denoise_post_REM(:, remove_cell_id) = [];
    all_denoise_post_REM = all_denoise_post_REM(:, s_idx);
    
    all_detrendRaw_pre_REM(:, remove_cell_id) = [];
    all_detrendRaw_pre_REM = all_detrendRaw_pre_REM(:, s_idx);
    
    all_detrendRaw_post_REM(:, remove_cell_id) = [];
    all_detrendRaw_post_REM = all_detrendRaw_post_REM(:, s_idx);
    
    % run the CSP extraction. 
    fprintf('Using Deconv/Denoise ("Grosmark et al. 2021")  data as an input\n')
    % pre
    REM_pre_data_in = []; 
    for ii = size(all_detrendRaw_pre_REM,2):-1:1
        m_ad(ii) = 1.5*mad([all_denoise_pre_REM(:,ii)- all_detrendRaw_pre_REM(:,ii)]);
        REM_pre_data_in(:,ii) = (all_deconv_pre_REM(:,ii)/m_ad(ii)) > 1.25;
    end
    REM_pre_data_in = REM_pre_data_in(:, s_idx); 
    
    % post
    REM_post_data_in = []; 
    for ii = size(all_detrendRaw_post_REM,2):-1:1
        m_ad(ii) = 1.5*mad([all_denoise_post_REM(:,ii)- all_detrendRaw_post_REM(:,ii)]);
        REM_post_data_in(:,ii) = (all_deconv_post_REM(:,ii)/m_ad(ii)) > 1.25;
    end
    REM_post_data_in = REM_post_data_in(:, s_idx); 
    
%     Csp =all_deconv_pre_REM./all_denoise_pre_REM;
%     Csp = Csp > 0.01;
%     REM_pre_data_in = Csp;
%     
%     % post
%     Csp =all_deconv_post_REM./all_denoise_post_REM;
%     Csp = Csp > 0.01;
%     REM_post_data_in = Csp;

    
else
    REM_pre_data_in = all_binary_pre_REM;
    REM_post_data_in = all_binary_post_REM; 
end

% sort the ms_trk_cut

f_list = fieldnames(ms_trk_cut); 

for ii = length(f_list):-1:1
    if size(ms_trk_cut.(f_list{ii}), 2) == length(S_neurons)
        disp(f_list{ii})
        ms_trk_cut.(f_list{ii}) = ms_trk_cut.(f_list{ii})(:, s_idx);
    elseif size(ms_trk_cut.(f_list{ii}), 3) == length(S_neurons)
        disp(['! ' f_list{ii}])
        ms_trk_cut.(f_list{ii}) = ms_trk_cut.(f_list{ii})(:,:, s_idx);
    end
    
end



%% load the place information
t_h5_dir = dir('tuning_*.h5');


PCs_properties = MS_h5_to_stuct(t_h5_dir(h5_idx).name); 




% load([info.subject '_' info.session '_PCs.mat'])


place.centroids = double(PCs_properties.peak_loc)'; 

% is it a place cell?
place.is = PCs_properties.p_value < 0.05;

place.map = PCs_properties.tuning_curves';

place.MI = PCs_properties.info;

place.peak_rate = PCs_properties.peak_val;



place.centroids(remove_cell_id)= [];
place.centroids =place.centroids(s_idx); 

place.is(remove_cell_id)= [];
place.is = place.is(s_idx); 


place.map(remove_cell_id,:)= [];
place.map =place.map(s_idx,:);

place.MI(remove_cell_id)= [];
place.MI = place.MI(s_idx); 

place.peak_rate(remove_cell_id)= [];
place.peak_rate = place.peak_rate(s_idx); 

% bin = 3; %80/size(place.map,2);
% p_bins = 0:bin:100;
% place.p_bins = p_bins(1:end)+bin/2;

place.bins = PCs_properties.bins; 
place.p_bins = PCs_properties.bins(1:end-1)+(mode(diff(PCs_properties.bins)))/2;

% see if there are any anxiety cells

% [~,p_sort] = sort(place.centroids);

%% get the initial assemblies

%     opts.threshold.method = 'binshuffling';
%     opts.Patterns.method = 'ICA';
%     opts.Patterns.number_of_iterations = 500;
%     opts.threshold.permutations_percentile= 95;
%     opts.threshold.number_of_permutations= 500;

A_temp = []; A_prog = []; wake_data = []; wake_tvec = [];
for iB = length(bin_s):-1:1
    
    [A_temp{iB}, A_proj{iB}, wake_data{iB}, wake_tvec{iB}, A_opts{iB}] = MS_PCA_ICA_only(ms_trk_cut, move_idx, bin_s(iB),method);
end


for iB = length(bin_s):-1:1
    fprintf('PCA-ICA detected %.0f assemblies using a %.2fs binsize\n', size(A_temp{iB},2), bin_s(iB))
end

%% if using the 'grosmark 2021' method, smooth the REM data as per the wake data
% if isfield(A_opts{1}, 'gauss_window')
%     
% 
%         gk = gausskernel(A_opts{1}.gauss_window,A_opts{1}.gauss_SD); 
%         gk = gk./A_opts{1}.binsize; % normalize by binsize
%         
%         REM_pre_data_in = conv2(REM_pre_data_in,gk,'same'); % convolve with gaussian window
%         
%         REM_post_data_in = conv2(REM_post_data_in,gk,'same'); % convolve with gaussian window
%         
% %         test = conv2(all_deconv_pre_REM(:,1)./all_denoise_pre_REM(:,1),gk,'same'); 
%         
%     
% end

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
for iB = length(bin_s):-1:1
    
    [P_loc{iB}] = MS_Asmbly_act_loc(P_proj{iB}, wake_tvec{iB}, behav, win_s, thresh, 2/bin_s(iB));
    
end


%% Load the REM data and compare it to the wake
cfg_ReAct.nShuff = 500;
cfg_ReAct.thresh = 99;

for iB = length(bin_s):-1:1
    % pre REM
    [REM_pre_proj{iB},REM_pre_stats{iB}, REM_pre_data{iB}, REM_pre_tvec{iB}, REM_Pre_shuff{iB}] = MS_Asmbly_ReAct(cfg_ReAct, REM_pre_data_in, P_temp{iB} ,ms_trk_cut,  bin_s(iB));
    
    % post REM
    [REM_post_proj{iB}, REM_post_stats{iB}, REM_post_data{iB}, REM_post_tvec{iB},REM_Post_shuff{iB}] = MS_Asmbly_ReAct(cfg_ReAct,REM_post_data_in, P_temp{iB} ,ms_trk_cut,  bin_s(iB));
    
end


%% get the cross correlation between assemblies
xc_bin = 1; 
t_max = 2; 

for iB = length(bin_s):-1:1
    
    [wake_zxcor{iB}, wake_zxcov{iB}, wake_xcor{iB},wake_xcov{iB}] = MS_Asmbly_xcor(P_proj{iB},wake_tvec{iB}, 5, xc_bin, t_max);
    
    [REM_pre_zxcor{iB}, REM_pre_zxcov{iB},REM_pre_xcor{iB},REM_pre_xcov{iB}] = MS_Asmbly_xcor(REM_pre_proj{iB},REM_pre_tvec{iB},REM_pre_stats{iB}.R_thresh , xc_bin, t_max);
    
    [REM_post_zxcor{iB}, REM_post_zxcov{iB},REM_post_xcor{iB}, REM_post_xcov{iB}] = MS_Asmbly_xcor(REM_post_proj{iB},REM_post_tvec{iB}, REM_post_stats{iB}.R_thresh, xc_bin, t_max);

    [REM_Pre_sig_CoOc{iB}, wake_sig_CoOc{iB}] = MS_Asmbly_CoAct_count(wake_zxcor{iB},REM_pre_zxcor{iB}, 1.96); 
    fprintf('Pre had %.0f (%.0f%%) significant assembly pairs using xcor of the %0.0f sig pairs found in wake \n', REM_Pre_sig_CoOc{iB}, (REM_Pre_sig_CoOc{iB} / wake_sig_CoOc{iB})*100, wake_sig_CoOc{iB})
    
    [REM_Post_sig_CoOc{iB}] = MS_Asmbly_CoAct_count(wake_zxcor{iB},REM_post_zxcor{iB}, 1.96); 
    fprintf('Post had %.0f (%.0f%%) significant assembliy pairs using xcor of the %0.0f sig pairs found in wake \n', REM_Post_sig_CoOc{iB}, (REM_Post_sig_CoOc{iB} / wake_sig_CoOc{iB})*100, wake_sig_CoOc{iB})

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
    info.move_thresh = move_thresh;
    info.method = method; 
    
    out{iB}.info = info;
    out{iB}.bins = bin_s(iB);
    % data
    out{iB}.behav = behav;
    out{iB}.move_idx = move_idx; 
    out{iB}.wake_data = wake_data{iB}; 
    out{iB}.wake_tvec = wake_tvec{iB}; 
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
    out{iB}.REM_Pre_cff = REM_pre_xcor{iB};
    out{iB}.REM_Pre_cffz = REM_pre_zxcor{iB};
    out{iB}.REM_Pre_nsig_cff = REM_Pre_sig_CoOc{iB}; 
    out{iB}.REM_Pre_psig_cff = (REM_Pre_sig_CoOc{iB}/wake_sig_CoOc{iB})*100;  

    
    out{iB}.REM_Post_proj = REM_post_proj{iB};
    out{iB}.REM_Post_stats = REM_post_stats{iB};
    out{iB}.REM_Post_data = REM_post_data{iB};
    out{iB}.REM_Post_tvec = REM_post_tvec{iB};
    out{iB}.REM_Post_shuff = REM_Post_shuff{iB};
    out{iB}.REM_Post_cff = REM_post_xcor{iB};
    out{iB}.REM_Post_cffz = REM_post_zxcor{iB};
    out{iB}.REM_Post_nsig_cff = REM_Post_sig_CoOc{iB}; 
    out{iB}.REM_Post_psig_cff = (REM_Post_sig_CoOc{iB}/wake_sig_CoOc{iB})*100;  
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
