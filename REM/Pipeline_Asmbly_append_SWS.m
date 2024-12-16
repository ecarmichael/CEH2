function out = Pipeline_Asmbly_append_SWS(fname, A_in)
%% Pipeline_Asmbly: provides a wrapper for running assembly and reactivation analyses using calcium data.


    % bin_s = A_in{1}.bins; 



rng(123, 'twister')


%% load the data

% load('ms_trk.mat')
% load('behav_DLC.mat')

load(fname, 'all_binary_pre_SW');
load(fname, 'all_binary_post_SW');

cell_id = 1:size(all_binary_pre_SW,2); 
remove_cell_id = find(~ismember(cell_id,double(A_in{1}.info.S_neurons)));

s_idx = []; 
for ii = length(A_in{1}.info.S_neurons):-1:1
   s_idx(ii) = find(A_in{1}.info.S_neurons(ii) == cell_id);
end

% trim the REM data as well
all_binary_pre_SW = all_binary_pre_SW(:, s_idx); 
all_binary_post_SW = all_binary_post_SW(:, s_idx); 



% % if using 'grosmark 2021' method, create the Csp 
% if strcmpi(method, 'grosmark')
%     % trim the REM data as well
%     all_deconv_pre_REM(:, remove_cell_id) = [];
%     all_deconv_pre_REM = all_deconv_pre_REM(:, s_idx);
%     
%     all_denoise_pre_REM(:, remove_cell_id) = [];
%     all_denoise_pre_REM = all_denoise_pre_REM(:, s_idx);
%     
%     % all_deconv_post_REM
%     
%     all_deconv_post_REM(:, remove_cell_id) = [];
%     all_deconv_post_REM = all_deconv_post_REM(:, s_idx);
%     
%     all_denoise_post_REM(:, remove_cell_id) = [];
%     all_denoise_post_REM = all_denoise_post_REM(:, s_idx);
%     
%     all_detrendRaw_pre_REM(:, remove_cell_id) = [];
%     all_detrendRaw_pre_REM = all_detrendRaw_pre_REM(:, s_idx);
%     
%     all_detrendRaw_post_REM(:, remove_cell_id) = [];
%     all_detrendRaw_post_REM = all_detrendRaw_post_REM(:, s_idx);
%     
%     % run the CSP extraction. 
%     fprintf('Using Deconv/Denoise ("Grosmark et al. 2021")  data as an input\n')
%     % pre
%     REM_pre_data_in = []; 
%     for ii = size(all_detrendRaw_pre_REM,2):-1:1
%         m_ad(ii) = 1.5*mad([all_denoise_pre_REM(:,ii)- all_detrendRaw_pre_REM(:,ii)]);
%         REM_pre_data_in(:,ii) = (all_deconv_pre_REM(:,ii)/m_ad(ii)) > 1.25;
%     end
%     REM_pre_data_in = REM_pre_data_in(:, s_idx); 
%     
%     % post
%     REM_post_data_in = []; 
%     for ii = size(all_detrendRaw_post_REM,2):-1:1
%         m_ad(ii) = 1.5*mad([all_denoise_post_REM(:,ii)- all_detrendRaw_post_REM(:,ii)]);
%         REM_post_data_in(:,ii) = (all_deconv_post_REM(:,ii)/m_ad(ii)) > 1.25;
%     end
%     REM_post_data_in = REM_post_data_in(:, s_idx); 
%     
% %     Csp =all_deconv_pre_REM./all_denoise_pre_REM;
% %     Csp = Csp > 0.01;
% %     REM_pre_data_in = Csp;
% %     
% %     % post
% %     Csp =all_deconv_post_REM./all_denoise_post_REM;
% %     Csp = Csp > 0.01;
% %     REM_post_data_in = Csp;
% 
%     
% else
    SW_pre_data_in = all_binary_pre_SW;
    SW_post_data_in = all_binary_post_SW; 
% end



%% place metrics

place = A_in{ii}.place; 



%% Load the REM data and compare it to the wake
cfg_ReAct.nShuff = 500;
cfg_ReAct.thresh = 99;
ms = []; 
cfg_ReAct.ms_fs = A_in{1}.info.ms_fs; 

for iB = length(A_in):-1:1
    % pre REM
    [SWS_pre_proj{iB},SWS_pre_stats{iB}, SWS_pre_data{iB}, SWS_pre_tvec{iB}, SWS_Pre_shuff{iB}] = MS_Asmbly_ReAct(cfg_ReAct, SW_pre_data_in, A_in{iB}.P_temp ,ms,  A_in{iB}.info.bin);
    
    % post REM
    [SWS_post_proj{iB}, SWS_post_stats{iB}, SWS_post_data{iB}, SWS_post_tvec{iB},SWS_Post_shuff{iB}] = MS_Asmbly_ReAct(cfg_ReAct, SW_post_data_in, A_in{iB}.P_temp ,ms,  A_in{iB}.info.bin);
    
end


%% get the cross correlation between assemblies
xc_bin = 1; 
t_max = 2; 

for iB = length(A_in):-1:1
    
    [wake_zxcor{iB}, wake_zxcov{iB}, wake_xcor{iB},wake_xcov{iB}] = MS_Asmbly_xcor(A_in{iB}.P_proj,A_in{iB}.wake_tvec, 5, xc_bin, t_max);
    
    [SWS_pre_zxcor{iB}, SWS_pre_zxcov{iB},SWS_pre_xcor{iB},SWS_pre_xcov{iB}] = MS_Asmbly_xcor(SWS_pre_proj{iB},SWS_pre_tvec{iB},SWS_pre_stats{iB}.R_thresh , xc_bin, t_max);
    
    [SWS_post_zxcor{iB}, SWS_post_zxcov{iB},SWS_post_xcor{iB}, SWS_post_xcov{iB}] = MS_Asmbly_xcor(SWS_post_proj{iB},SWS_post_tvec{iB}, SWS_post_stats{iB}.R_thresh, xc_bin, t_max);

    [SWS_Pre_sig_CoOc{iB}, wake_sig_CoOc{iB}] = MS_Asmbly_CoAct_count(wake_zxcor{iB},SWS_pre_zxcor{iB}, 1.96); 
    fprintf('Pre had %.0f (%.0f%%) significant assembly pairs using xcor of the %0.0f sig pairs found in wake \n', SWS_Pre_sig_CoOc{iB}, (SWS_Pre_sig_CoOc{iB} / wake_sig_CoOc{iB})*100, wake_sig_CoOc{iB})
    
    [SWS_Post_sig_CoOc{iB}] = MS_Asmbly_CoAct_count(wake_zxcor{iB},SWS_post_zxcor{iB}, 1.96); 
    fprintf('Post had %.0f (%.0f%%) significant assembliy pairs using xcor of the %0.0f sig pairs found in wake \n', SWS_Post_sig_CoOc{iB}, (SWS_Post_sig_CoOc{iB} / wake_sig_CoOc{iB})*100, wake_sig_CoOc{iB})

end


%% reactivation strength

for iB = length(A_in):-1:1
    
    for ii = size(A_in{iB}.P_proj,1):-1:1
        
        SWS_ReAct{iB}(ii) = mean(SWS_post_proj{iB}(ii,:)) - mean(SWS_pre_proj{iB}(ii,:));
    end
end

%% collect the outputs
for iB = length(A_in):-1:1
    out{iB} = A_in{iB}; 

    
    out{iB}.SWS_pre_in = SW_pre_data_in;
    out{iB}.SWS_post_in = SW_post_data_in;
    
    % reactivations
    out{iB}.SWS_ReAct = SWS_ReAct{iB};
    
    out{iB}.SWS_Pre_proj = SWS_pre_proj{iB};
    out{iB}.SWS_Pre_stats = SWS_pre_stats{iB};
    out{iB}.SWS_Pre_data = SWS_pre_data{iB};
    out{iB}.SWS_Pre_tvec = SWS_pre_tvec{iB};
    out{iB}.SWS_Pre_shuff = SWS_Pre_shuff{iB};
    out{iB}.SWS_Pre_cff = SWS_pre_xcor{iB};
    out{iB}.SWS_Pre_cffz = SWS_pre_zxcor{iB};
    out{iB}.SWS_Pre_nsig_cff = SWS_Pre_sig_CoOc{iB}; 
    out{iB}.SWS_Pre_psig_cff = (SWS_Pre_sig_CoOc{iB}/wake_sig_CoOc{iB})*100;  

    
    out{iB}.SWS_Post_proj = SWS_post_proj{iB};
    out{iB}.SWS_Post_stats = SWS_post_stats{iB};
    out{iB}.SWS_Post_data = SWS_post_data{iB};
    out{iB}.SWS_Post_tvec = SWS_post_tvec{iB};
    out{iB}.SWS_Post_shuff = SWS_Post_shuff{iB};
    out{iB}.SWS_Post_cff = SWS_post_xcor{iB};
    out{iB}.SWS_Post_cffz = SWS_post_zxcor{iB};
    out{iB}.SWS_Post_nsig_cff = SWS_Post_sig_CoOc{iB}; 
    out{iB}.SWS_Post_psig_cff = (SWS_Post_sig_CoOc{iB}/wake_sig_CoOc{iB})*100;  
end

