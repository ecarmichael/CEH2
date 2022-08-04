%% JC pairwise correlation sandbox

data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_18_2021_pv1252_HATD1';

cd(data_dir)

sess_list = {'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_18_2021_pv1252_HATD1'}; 

parts = strsplit(sess_list{1}, filesep); 
sub = parts{end-1};
sub = [lower(sub(1:2)) sub(3:end)]; 
sess_parts = strsplit(parts{end}, '_'); 
sess = sess_parts{end}; 

% load the selective cell classificiation
selective = load([strjoin(parts(1:4), filesep) filesep '2.Selective_cell' filesep sub '_' sess '.mat']); 
selective = selective.ans; 
%% load some data
load('ms_trk.mat')

trk_detrendraw = ms_trk.detrendRaw; 

load('all_detrendRaw_pre_REM.mat')
load('all_detrendRaw_post_REM.mat')

clear ms


%% compute the pairwise correlation across each segment

% pre REM
[pre_rho, pre_pval] = corr(all_detrendRaw_pre_REM, 'rows', 'pairwise'); 
% trk
[trk_rho, trk_pval] = corr(trk_detrendraw, 'rows', 'pairwise'); 
% post REM
[post_rho, post_pval] = corr(all_detrendRaw_post_REM, 'rows', 'pairwise'); 

% convert to dissimilarity 
pre_rho_corr = pre_rho - eye(size(pre_rho)); 
pre_rho_dis = 1 - pre_rho_corr(find(pre_rho_corr))';

trk_rho_corr = trk_rho - eye(size(trk_rho)); 
trk_rho_dis = 1 - trk_rho_corr(find(trk_rho_corr))'; 

post_rho_corr = post_rho - eye(size(post_rho)); 
post_rho_dis = 1 - post_rho_corr(find(post_rho_corr))';  

%% count the significantly correlated pairs
idx_mat = triu(ones(size(pre_pval)), 1); % index for upper half of matrix, above diagonal. 

pre_sig_cells = NaN(size(pre_pval));
trk_sig_cells = NaN(size(trk_pval));
post_sig_cells = NaN(size(post_pval));

for ii = 1:size(idx_mat,1)
    for jj = 1:size(idx_mat,2)
        if idx_mat(ii, jj)
            pre_sig_cells(ii, jj) = pre_pval(ii, jj);
            trk_sig_cells(ii, jj) = trk_pval(ii, jj);
            post_sig_cells(ii, jj) = post_pval(ii, jj);
        end
    end
end

pre_sig_cells(pre_sig_cells > 0.05) = NaN; 
trk_sig_cells(trk_sig_cells > 0.05) = NaN; 
post_sig_cells(post_sig_cells > 0.05) = NaN; 

pre_nsig =(sum(~isnan(pre_sig_cells), 'all', 'omitnan')/ numel(pre_sig_cells))*100; 
trk_nsig = (sum(~isnan(trk_sig_cells), 'all', 'omitnan')/ numel(trk_sig_cells))*100; 
post_nsig = (sum(~isnan(post_sig_cells), 'all', 'omitnan')/ numel(post_sig_cells))*100;

fprintf('Pre sig: %2.2f%% \n', pre_nsig)
fprintf('trk sig: %2.2f%% \n', trk_nsig)
fprintf('Post sig: %2.2f%% \n',  post_nsig)

%% plot 

figure(1)
clf
% corr mat
subplot(2,3,1)
imagesc(pre_rho_corr); 
set(gca, 'XAxisLocation', 'top')
title('Pre REM')
% colorbar; 
caxis([0 1]); 

subplot(2,3,2)
imagesc(trk_rho); 
set(gca, 'XAxisLocation', 'top')
title('Track')
% colorbar; 
caxis([0 1]); 

subplot(2,3,3)
imagesc(post_rho); 
set(gca, 'XAxisLocation', 'top')
title('Post REM')
% colorbar; 
caxis([0 1]); 

% pval mat
subplot(2,3,4)
imagesc(pre_pval); 
set(gca, 'XAxisLocation', 'top')
title(['Pre REM p val (' num2str(pre_nsig, '%2.2f') '%)'])
% colorbar; 
caxis([0 .05]); 

subplot(2,3,5)
imagesc(trk_pval); 
set(gca, 'XAxisLocation', 'top')
title(['Track p val (' num2str(trk_nsig, '%2.2f') '%)'])
% colorbar; 
caxis([0 .05]); 

subplot(2,3,6)
imagesc(post_pval); 
set(gca, 'XAxisLocation', 'top')
title(['Post REM p val (' num2str(post_nsig, '%2.2f') '%)'])
% colorbar; 
caxis([0 .05]); 
% 
% % dis mat
% subplot(3,3,7)
% imagesc(pre_rho_dis); 
% set(gca, 'XAxisLocation', 'top')
% title('Pre REM Dis')
% % colorbar; 
% caxis([0 2]); 
% 
% subplot(3,3,8)
% imagesc(trk_rho_dis); 
% set(gca, 'XAxisLocation', 'top')
% title('Track Dis')
% % colorbar; 
% caxis([0 2]); 
% 
% subplot(3,3,9)
% imagesc(post_rho_dis); 
% set(gca, 'XAxisLocation', 'top')
% title('Post REM Dis')
% % colorbar; 
% caxis([0 2]); 


%% seletive cells only


if exist('selective')
    all_detrendRaw_pre_REM_L = all_detrendRaw_pre_REM(:, ~isnan(selective.Left_Selective)); 
    all_detrendRaw_pre_REM_R = all_detrendRaw_pre_REM(:, ~isnan(selective.Right_Selective)); 
    
    trk_detrendraw_L = trk_detrendraw(:, ~isnan(selective.Left_Selective)); 
    trk_detrendraw_R = trk_detrendraw(:, ~isnan(selective.Right_Selective)); 
    
    all_detrendRaw_post_REM_L = all_detrendRaw_post_REM(:, ~isnan(selective.Left_Selective)); 
    all_detrendRaw_post_REM_R = all_detrendRaw_post_REM(:, ~isnan(selective.Right_Selective)); 
    
    
    % pre REM
[pre_rho_L, pre_pval_L] = corr(all_detrendRaw_pre_REM_L, 'rows', 'pairwise'); 
[pre_rho_R, pre_pval_R] = corr(all_detrendRaw_pre_REM_R, 'rows', 'pairwise'); 

% trk
[trk_rho_L, trk_pval_L] = corr(trk_detrendraw_L, 'rows', 'pairwise'); 
[trk_rho_R, trk_pval_R] = corr(trk_detrendraw_R, 'rows', 'pairwise'); 

% post REM
[post_rho_L, post_pval_L] = corr(all_detrendRaw_post_REM_L, 'rows', 'pairwise'); 
[post_rho_R, post_pval_R] = corr(all_detrendRaw_post_REM_R, 'rows', 'pairwise'); 


% get sig for L selective. 
idx_mat_L = triu(ones(size(pre_pval_L)), 1); % index for upper half of matrix, above diagonal. 

pre_L_sig_cells = NaN(size(idx_mat_L));
trk_L_sig_cells = NaN(size(idx_mat_L));
post_L_sig_cells = NaN(size(idx_mat_L));

for ii = 1:size(idx_mat_L,1)
    for jj = 1:size(idx_mat_L,2)
        if idx_mat_L(ii, jj)
            pre_L_sig_cells(ii, jj) = pre_pval_L(ii, jj);
            trk_L_sig_cells(ii, jj) = trk_pval_L(ii, jj);
            post_L_sig_cells(ii, jj) = post_pval_L(ii, jj);
        end
    end
end

pre_L_sig_cells(pre_L_sig_cells > 0.05) = NaN; 
trk_L_sig_cells(trk_L_sig_cells > 0.05) = NaN; 
post_L_sig_cells(post_L_sig_cells > 0.05) = NaN; 

pre_L_nsig =(sum(~isnan(pre_L_sig_cells), 'all', 'omitnan')/ numel(pre_L_sig_cells))*100; 
trk_L_nsig = (sum(~isnan(trk_L_sig_cells), 'all', 'omitnan')/ numel(trk_L_sig_cells))*100; 
post_L_nsig = (sum(~isnan(post_L_sig_cells), 'all', 'omitnan')/ numel(post_L_sig_cells))*100;

fprintf('Pre sig: %2.2f%% \n', pre_L_nsig)
fprintf('trk sig: %2.2f%% \n', trk_L_nsig)
fprintf('Post sig: %2.2f%% \n',  post_L_nsig)

% right seletive
idx_mat_R = triu(ones(size(pre_pval_R)), 1); % index for upper half of matrix, above diagonal. 

pre_R_sig_cells = NaN(size(idx_mat_R));
trk_R_sig_cells = NaN(size(idx_mat_R));
post_R_sig_cells = NaN(size(idx_mat_R));

for ii = 1:size(idx_mat_R,1)
    for jj = 1:size(idx_mat_R,2)
        if idx_mat_R(ii, jj)
            pre_R_sig_cells(ii, jj) = pre_pval_R(ii, jj);
            trk_R_sig_cells(ii, jj) = trk_pval_R(ii, jj);
            post_R_sig_cells(ii, jj) = post_pval_R(ii, jj);
        end
    end
end

pre_R_sig_cells(pre_R_sig_cells > 0.05) = NaN; 
trk_R_sig_cells(trk_R_sig_cells > 0.05) = NaN; 
post_R_sig_cells(post_R_sig_cells > 0.05) = NaN; 

pre_R_nsig =(sum(~isnan(pre_R_sig_cells), 'all', 'omitnan')/ numel(pre_R_sig_cells))*100; 
trk_R_nsig = (sum(~isnan(trk_R_sig_cells), 'all', 'omitnan')/ numel(trk_R_sig_cells))*100; 
post_R_nsig = (sum(~isnan(post_R_sig_cells), 'all', 'omitnan')/ numel(post_R_sig_cells))*100;

fprintf('Pre sig: %2.2f%% \n', pre_R_nsig)
fprintf('trk sig: %2.2f%% \n', trk_R_nsig)
fprintf('Post sig: %2.2f%% \n',  post_R_nsig)
    
    
% plot
figure(2)
clf
% corr mat
subplot(4,3,1)
imagesc(pre_rho_L); 
set(gca, 'XAxisLocation', 'top')
title('Pre REM L selective')
% colorbar; 
caxis([0 1]); 

subplot(4,3,2)
imagesc(trk_rho_L); 
set(gca, 'XAxisLocation', 'top')
title('Track L')
% colorbar; 
caxis([0 1]); 

subplot(4,3,3)
imagesc(post_rho_L); 
set(gca, 'XAxisLocation', 'top')
title('Post REM L')
% colorbar; 
caxis([0 1]); 

% pval mat
subplot(4,3,4)
imagesc(pre_pval_L); 
set(gca, 'XAxisLocation', 'top')
title(['Pre REM L p val (' num2str(pre_L_nsig, '%2.2f') '%)'])
% colorbar; 
caxis([0 .05]); 

subplot(4,3,5)
imagesc(trk_pval_L); 
set(gca, 'XAxisLocation', 'top')
title(['Track L p val (' num2str(trk_L_nsig, '%2.2f') '%)'])
% colorbar; 
caxis([0 .05]); 

subplot(4,3,6)
imagesc(post_pval_L); 
set(gca, 'XAxisLocation', 'top')
title(['Post REM L p val (' num2str(post_L_nsig, '%2.2f') '%)'])
% colorbar; 
caxis([0 .05]); 

% right
subplot(4,3,7)
imagesc(pre_rho_R); 
set(gca, 'XAxisLocation', 'top')
title('Pre REM R selective')
% colorbar; 
caxis([0 1]); 

subplot(4,3,8)
imagesc(trk_rho_R); 
set(gca, 'XAxisLocation', 'top')
title('Track R')
% colorbar; 
caxis([0 1]); 

subplot(4,3,9)
imagesc(post_rho_R); 
set(gca, 'XAxisLocation', 'top')
title('Post REM R')
% colorbar; 
caxis([0 1]); 

% pval mat
subplot(4,3,10)
imagesc(pre_pval_R); 
set(gca, 'XAxisLocation', 'top')
title(['Pre REM R p val (' num2str(pre_R_nsig, '%2.2f') '%)'])
% colorbar; 
caxis([0 .05]); 

subplot(4,3,11)
imagesc(trk_pval_R); 
set(gca, 'XAxisLocation', 'top')
title(['Track R p val (' num2str(trk_R_nsig, '%2.2f') '%)'])
% colorbar; 
caxis([0 .05]); 

subplot(4,3,12)
imagesc(post_pval_R); 
set(gca, 'XAxisLocation', 'top')
title(['Post REM R p val (' num2str(post_R_nsig, '%2.2f') '%)'])
% colorbar; 
caxis([0 .05]); 

    
    
end
