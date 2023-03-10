%% JC pairwise correlation sandbox

inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Correlation_plots_EC';

data_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_18_2021_pv1252_HATD1';

cd(data_dir)

%%
sess_list = {'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1060\7_15_2019_PV1060_LTD1';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1060\7_19_2019_PV1060_LTD5';...
    'E:\Jisoo_Project\Inter\PV1060\11_19_2019_PV1060_HATD1';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1060\11_26_2019_PV1060_HATSwitch';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_8_2019_PV1069_LTD1';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\7_12_2019_PV1069_LTD5';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_18_2019_PV1069_HATD5';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1069\10_22_2019_PV1069_HATSwitch';...
    'E:\Jisoo_Project\Inter\PV1191\5_19_2021_PV1191_HATD1';...
    'E:\Jisoo_Project\Inter\PV1191\5_23_2021_PV1191_HATD5';...
    'E:\Jisoo_Project\Inter\PV1191\5_25_2021_PV1191_HATDS';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1192\4_17_2021_PV1192_HATD1';...
    'E:\Jisoo_Project\Inter\PV1192\4_21_2021_PV1192_HATD5';...
    'E:\Jisoo_Project\Inter\PV1192\4_23_2021_PV1192_HATDSwitch';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_16_2021_pv1252_LTD5'; ...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_18_2021_pv1252_HATD1';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_22_2021_pv1252_HATD5';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1252\11_24_2021_pv1252_HATDSwitch';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1254\11_13_2021_pv1254_LTD1';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1254\11_17_2021_pv1254_LTD5';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1254\11_19_2021_pv1254_HATD1';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1254\11_23_2021_pv1254_HATD5';...
    'C:\Users\ecarm\Dropbox (Williams Lab)\JisooProject2020\2020_Results_aftercutting\Across_episodes\Inter\PV1254\11_25_2021_pv1254_HATDSwitch';...
    };

for iSess = 1:length(sess_list)
    cd(sess_list{iSess});
    parts = strsplit(sess_list{iSess}, filesep);
    sub = parts{end-1};
    sub = [lower(sub(1:2)) sub(3:end)];
    sess_parts = strsplit(parts{end}, '_');
    sess = sess_parts{end};
    
    if ~exist(['C:\Users\ecarm\Dropbox (Williams Lab)\' filesep '2.Selective_cell' filesep 'Jimenez_methods' filesep sub filesep sess filesep 'Out.mat']) || ~exist(['C:\Users\ecarm\Dropbox (Williams Lab)\4.PlaceCell'  filesep sub filesep sess filesep 'spatial_analysis_classif.mat'])
        continue
    end
    % load the selective cell classificiation
    selective = load(['C:\Users\ecarm\Dropbox (Williams Lab)\' filesep '2.Selective_cell' filesep 'Jimenez_methods' filesep sub filesep sess filesep 'Out.mat']);
    open = selective.ans.Jime_OpenCell_SD1';
    closed = selective.ans.Jime_ClosedCell_SD1';
    
    %load place selective
    place_data = load(['C:\Users\ecarm\Dropbox (Williams Lab)\4.PlaceCell'  filesep sub filesep sess filesep 'spatial_analysis_classif.mat']);
    place = place_data.SA.WholePlaceCell;
    
    for iP = length(place):-1:1
        if iscell(place_data.SA.PlaceFieldCentroid{place(iP),1})
            place_centeroid(iP) = place_data.SA.PlaceFieldCentroid{place(iP),1}{1,1}(1);
        elseif iscell(place_data.SA.PlaceFieldCentroid{place(iP),2})
            place_centeroid(iP) = place_data.SA.PlaceFieldCentroid{place(iP),2}{1,1}(1);
        elseif iscell(place_data.SA.PlaceFieldCentroid{place(iP),3})
            place_centeroid(iP) = place_data.SA.PlaceFieldCentroid{place(iP),3}{1,1}(1);
        end
    end
    
    % sort cells based on pl centroid
    [cent_s, s_idx] = sort(place_centeroid);
    place_s = place(s_idx);
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
    maximize
    
    % corr mat
    subplot(2,3,1)
    imagesc(pre_rho_corr);
    set(gca, 'XAxisLocation', 'top')
    title('Pre REM')
    % colorbar;
    caxis([0 1]);
    
    subplot(2,3,2)
    imagesc(trk_rho_corr);
    set(gca, 'XAxisLocation', 'top')
    title('Track')
    % colorbar;
    caxis([0 1]);
    
    subplot(2,3,3)
    imagesc(post_rho_corr);
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
    
    % save figures
    saveas(gcf, [inter_dir filesep sub '_' sess '_all_corr.fig'])
    saveas(gcf, [inter_dir filesep sub '_' sess '_all_corr.png'])
    
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
    
    % collect the data
    all_out.(sub).(sess).all.pre_nsig = pre_nsig;
    all_out.(sub).(sess).all.trk_nsig = trk_nsig;
    all_out.(sub).(sess).all.post_nsig = post_nsig;
    all_out.(sub).(sess).all.pre_sig_ind = ~isnan(pre_sig_cells);
    all_out.(sub).(sess).all.trk_sig_ind = ~isnan(trk_sig_cells);
    all_out.(sub).(sess).all.post_sig_ind = ~isnan(post_sig_cells);
    
    %% seletive cells only
    
    
    if exist('open')
        
        % pre REM
        [pre_rho_O, pre_pval_O] = corr(all_detrendRaw_pre_REM(:, open), 'rows', 'pairwise');
        [pre_rho_C, pre_pval_C] = corr(all_detrendRaw_pre_REM(:, closed), 'rows', 'pairwise');
        
        % trk
        [trk_rho_O, trk_pval_O] = corr(trk_detrendraw(:, open), 'rows', 'pairwise');
        [trk_rho_C, trk_pval_C] = corr(trk_detrendraw(:, closed), 'rows', 'pairwise');
        
        % post REM
        [post_rho_O, post_pval_O] = corr(all_detrendRaw_post_REM(:, open), 'rows', 'pairwise');
        [post_rho_C, post_pval_C] = corr(all_detrendRaw_post_REM(:, closed), 'rows', 'pairwise');
        
        
        % get sig for L selective.
        idx_mat_O = triu(ones(size(pre_pval_O)), 1); % index for upper half of matrix, above diagonal.
        
        pre_O_sig_cells = NaN(size(idx_mat_O));
        trk_O_sig_cells = NaN(size(idx_mat_O));
        post_O_sig_cells = NaN(size(idx_mat_O));
        
        for ii = 1:size(idx_mat_O,1)
            for jj = 1:size(idx_mat_O,2)
                if idx_mat_O(ii, jj)
                    pre_O_sig_cells(ii, jj) = pre_pval_O(ii, jj);
                    trk_O_sig_cells(ii, jj) = trk_pval_O(ii, jj);
                    post_O_sig_cells(ii, jj) = post_pval_O(ii, jj);
                end
            end
        end
        
        pre_O_sig_cells(pre_O_sig_cells > 0.05) = NaN;
        trk_O_sig_cells(trk_O_sig_cells > 0.05) = NaN;
        post_O_sig_cells(post_O_sig_cells > 0.05) = NaN;
        
        pre_O_nsig =(sum(~isnan(pre_O_sig_cells), 'all', 'omitnan')/ numel(pre_O_sig_cells))*100;
        trk_O_nsig = (sum(~isnan(trk_O_sig_cells), 'all', 'omitnan')/ numel(trk_O_sig_cells))*100;
        post_O_nsig = (sum(~isnan(post_O_sig_cells), 'all', 'omitnan')/ numel(post_O_sig_cells))*100;
        
        fprintf('Pre Open sig: %2.2f%% \n', pre_O_nsig)
        fprintf('trk Open sig: %2.2f%% \n', trk_O_nsig)
        fprintf('Post Open sig: %2.2f%% \n',  post_O_nsig)
        
        % right seletive
        idx_mat_C = triu(ones(size(pre_pval_C)), 1); % index for upper half of matrix, above diagonal.
        
        pre_C_sig_cells = NaN(size(idx_mat_C));
        trk_C_sig_cells = NaN(size(idx_mat_C));
        post_C_sig_cells = NaN(size(idx_mat_C));
        
        for ii = 1:size(idx_mat_C,1)
            for jj = 1:size(idx_mat_C,2)
                if idx_mat_C(ii, jj)
                    pre_C_sig_cells(ii, jj) = pre_pval_C(ii, jj);
                    trk_C_sig_cells(ii, jj) = trk_pval_C(ii, jj);
                    post_C_sig_cells(ii, jj) = post_pval_C(ii, jj);
                end
            end
        end
        
        pre_C_sig_cells(pre_C_sig_cells > 0.05) = NaN;
        trk_C_sig_cells(trk_C_sig_cells > 0.05) = NaN;
        post_C_sig_cells(post_C_sig_cells > 0.05) = NaN;
        
        pre_C_nsig =(sum(~isnan(pre_C_sig_cells), 'all', 'omitnan')/ numel(pre_C_sig_cells))*100;
        trk_C_nsig = (sum(~isnan(trk_C_sig_cells), 'all', 'omitnan')/ numel(trk_C_sig_cells))*100;
        post_C_nsig = (sum(~isnan(post_C_sig_cells), 'all', 'omitnan')/ numel(post_C_sig_cells))*100;
        
        fprintf('Pre Closed sig: %2.2f%% \n', pre_C_nsig)
        fprintf('trk Closed sig: %2.2f%% \n', trk_C_nsig)
        fprintf('Post Closed sig: %2.2f%% \n',  post_C_nsig)
        
        
        % plot
        figure(2)
        clf
        maximize
        
        % corr mat
        subplot(4,3,1)
        imagesc(pre_rho_O);
        set(gca, 'XAxisLocation', 'top')
        title('Pre REM Open selective')
        % colorbar;
        caxis([0 1]);
        
        subplot(4,3,2)
        imagesc(trk_rho_O);
        set(gca, 'XAxisLocation', 'top')
        title('Track \color{yellow}Open')
        % colorbar;
        caxis([0 1]);
        
        subplot(4,3,3)
        imagesc(post_rho_O);
        set(gca, 'XAxisLocation', 'top')
        title('Post REM Open')
        % colorbar;
        caxis([0 1]);
        
        % pval mat
        subplot(4,3,4)
        imagesc(pre_pval_O);
        set(gca, 'XAxisLocation', 'top')
        title(['Pre REM Open p val (' num2str(pre_O_nsig, '%2.2f') '%)'])
        % colorbar;
        caxis([0 .05]);
        
        subplot(4,3,5)
        imagesc(trk_pval_O);
        set(gca, 'XAxisLocation', 'top')
        title(['Track Open p val (' num2str(trk_O_nsig, '%2.2f') '%)'])
        % colorbar;
        caxis([0 .05]);
        
        subplot(4,3,6)
        imagesc(post_pval_O);
        set(gca, 'XAxisLocation', 'top')
        title(['Post REM Open p val (' num2str(post_O_nsig, '%2.2f') '%)'])
        % colorbar;
        caxis([0 .05]);
        
        % right
        subplot(4,3,7)
        imagesc(pre_rho_C);
        set(gca, 'XAxisLocation', 'top')
        title('Pre REM Closed selective')
        % colorbar;
        caxis([0 1]);
        
        subplot(4,3,8)
        imagesc(trk_rho_C);
        set(gca, 'XAxisLocation', 'top')
        title('Track \color{blue}Closed')
        % colorbar;
        caxis([0 1]);
        
        subplot(4,3,9)
        imagesc(post_rho_C);
        set(gca, 'XAxisLocation', 'top')
        title('Post REM Closed')
        % colorbar;
        caxis([0 1]);
        
        % pval mat
        subplot(4,3,10)
        imagesc(pre_pval_C);
        set(gca, 'XAxisLocation', 'top')
        title(['Pre REM Closed p val (' num2str(pre_C_nsig, '%2.2f') '%)'])
        % colorbar;
        caxis([0 .05]);
        
        subplot(4,3,11)
        imagesc(trk_pval_C);
        set(gca, 'XAxisLocation', 'top')
        title(['Track Closed p val (' num2str(trk_C_nsig, '%2.2f') '%)'])
        % colorbar;
        caxis([0 .05]);
        
        subplot(4,3,12)
        imagesc(post_pval_C);
        set(gca, 'XAxisLocation', 'top')
        title(['Post REM Closed p val (' num2str(post_C_nsig, '%2.2f') '%)'])
        % colorbar;
        caxis([0 .05]);
        
        % save figures
        saveas(gcf, [inter_dir filesep sub '_' sess '_OC_selective_corr.fig'])
        saveas(gcf, [inter_dir filesep sub '_' sess '_OC_selective_corr.png'])
        
        % collect the data
        all_out.(sub).(sess).open.pre_nsig = pre_O_nsig;
        all_out.(sub).(sess).open.trk_nsig = trk_O_nsig;
        all_out.(sub).(sess).open.post_nsig = post_O_nsig;
        all_out.(sub).(sess).open.pre_sig_ind = ~isnan(pre_O_sig_cells);
        all_out.(sub).(sess).open.trk_sig_ind = ~isnan(trk_O_sig_cells);
        all_out.(sub).(sess).open.post_sig_ind = ~isnan(post_O_sig_cells);
        
        all_out.(sub).(sess).close.pre_nsig = pre_C_nsig;
        all_out.(sub).(sess).close.trk_nsig = trk_C_nsig;
        all_out.(sub).(sess).close.post_nsig = post_C_nsig;
        all_out.(sub).(sess).close.pre_sig_ind = ~isnan(pre_C_sig_cells);
        all_out.(sub).(sess).close.trk_sig_ind = ~isnan(trk_C_sig_cells);
        all_out.(sub).(sess).close.post_sig_ind = ~isnan(post_C_sig_cells);
        
    end
    
    
    %% place only
    
    if exist('place_s')
        
        mid_9p5_trk_idx = nearest_idx(9.5, cent_s);
        mid_10p5_trk_idx = nearest_idx(10.5, cent_s);
        
        
        % pre REM
        [pre_rho, pre_pval] = corr(all_detrendRaw_pre_REM(:, place), 'rows', 'pairwise');
        % trk
        [trk_rho, trk_pval] = corr(trk_detrendraw(:, place), 'rows', 'pairwise');
        % post REM
        [post_rho, post_pval] = corr(all_detrendRaw_post_REM(:, place), 'rows', 'pairwise');
        
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
        
        figure(3)
        clf
        maximize
        
        % corr mat
        subplot(2,3,1)
        imagesc(pre_rho_corr);
        set(gca, 'XAxisLocation', 'top')
        title('Pre REM')
        % colorbar;
        caxis([0 1]);
        yline(mid_9p5_trk_idx,'-',  '<9.5cm', 'color', [0.9153    0.2816    0.2878]);
        yline(mid_10p5_trk_idx,'-',  '>10.5cm', 'color', [0.9153    0.2816    0.2878], 'LabelVerticalAlignment', 'bottom');
        
        subplot(2,3,2)
        imagesc(trk_rho_corr);
        set(gca, 'XAxisLocation', 'top')
        title('Track \color{red}Place')
        % colorbar;
        caxis([0 1]);
        yline(mid_9p5_trk_idx,'-',  '<9.5cm', 'color', [0.9153    0.2816    0.2878]);
        yline(mid_10p5_trk_idx,'-',  '>10.5cm', 'color', [0.9153    0.2816    0.2878], 'LabelVerticalAlignment', 'bottom');
        
        
        subplot(2,3,3)
        imagesc(post_rho_corr);
        set(gca, 'XAxisLocation', 'top')
        title('Post REM')
        % colorbar;
        caxis([0 1]);
        yline(mid_9p5_trk_idx,'-',  '<9.5cm', 'color', [0.9153    0.2816    0.2878]);
        yline(mid_10p5_trk_idx,'-',  '>10.5cm', 'color', [0.9153    0.2816    0.2878], 'LabelVerticalAlignment', 'bottom');
        
        % pval mat
        subplot(2,3,4)
        imagesc(pre_pval);
        set(gca, 'XAxisLocation', 'top')
        title(['Pre REM p val (' num2str(pre_nsig, '%2.2f') '%)'])
        % colorbar;
        caxis([0 .05]);
        yline(mid_9p5_trk_idx,'-',  '<9.5cm', 'color', [0.9153    0.2816    0.2878]);
        yline(mid_10p5_trk_idx,'-',  '>10.5cm', 'color', [0.9153    0.2816    0.2878], 'LabelVerticalAlignment', 'bottom');
        
        
        subplot(2,3,5)
        imagesc(trk_pval);
        set(gca, 'XAxisLocation', 'top')
        title(['Track p val (' num2str(trk_nsig, '%2.2f') '%)'])
        % colorbar;
        caxis([0 .05]);
        yline(mid_9p5_trk_idx,'-',  '<9.5cm', 'color', [0.9153    0.2816    0.2878]);
        yline(mid_10p5_trk_idx,'-',  '>10.5cm', 'color', [0.9153    0.2816    0.2878], 'LabelVerticalAlignment', 'bottom');
        
        subplot(2,3,6)
        imagesc(post_pval);
        set(gca, 'XAxisLocation', 'top')
        title(['Post REM p val (' num2str(post_nsig, '%2.2f') '%)'])
        % colorbar;
        caxis([0 .05]);
        yline(mid_9p5_trk_idx,'-',  '<9.5cm', 'color', [0.9153    0.2816    0.2878]);
        yline(mid_10p5_trk_idx,'-',  '>10.5cm', 'color', [0.9153    0.2816    0.2878], 'LabelVerticalAlignment', 'bottom');
        
        % save figures
        saveas(gcf, [inter_dir filesep sub '_' sess '_place_corr.fig'])
        saveas(gcf, [inter_dir filesep sub '_' sess '_place_corr.png'])
        
        
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
        
        % collect the data
        all_out.(sub).(sess).place.pre_nsig = pre_nsig;
        all_out.(sub).(sess).place.trk_nsig = trk_nsig;
        all_out.(sub).(sess).place.post_nsig = post_nsig;
        all_out.(sub).(sess).place.pre_sig_ind = ~isnan(pre_sig_cells);
        all_out.(sub).(sess).place.trk_sig_ind = ~isnan(trk_sig_cells);
        all_out.(sub).(sess).place.post_sig_ind = ~isnan(post_sig_cells);
        
    end
    
    close all
    
    save([inter_dir filesep 'All_corr_out.mat'], 'all_out')
    clearvars -except sess_list inter_dir iSess all_out
end% loop over sessions

%% collect data and make a summary plot
subs = fieldnames(all_out);

all.LTD1_pre = NaN(size(subs)); all.LTD1_trk = NaN(size(subs)); all.LTD1_post = NaN(size(subs));
all.LTD5_pre = NaN(size(subs)); all.LTD5_trk = NaN(size(subs)); all.LTD5_post = NaN(size(subs));
all.HATD1_pre = NaN(size(subs)); all.HATD1_trk = NaN(size(subs)); all.HATD1_post = NaN(size(subs));
all.HATD5_pre = NaN(size(subs)); all.HATD5_trk = NaN(size(subs)); all.HATD5_post = NaN(size(subs));
all.HATDS_pre = NaN(size(subs)); all.HATDS_trk = NaN(size(subs)); all.HATDS_post = NaN(size(subs));
open = all; closed = all; place = all;

for ii = 1:length(subs)
    if isfield(all_out.(subs{ii}), 'LTD1')
        all.LTD1_pre(ii) = all_out.(subs{ii}).LTD1.all.pre_nsig;
        open.LTD1_pre(ii) = all_out.(subs{ii}).LTD1.open.pre_nsig;
        closed.LTD1_pre(ii) = all_out.(subs{ii}).LTD1.close.pre_nsig;
        place.LTD1_pre(ii) = all_out.(subs{ii}).LTD1.place.pre_nsig;
        
        all.LTD1_trk(ii) = all_out.(subs{ii}).LTD1.all.trk_nsig;
        open.LTD1_trk(ii) = all_out.(subs{ii}).LTD1.open.trk_nsig;
        closed.LTD1_trk(ii) = all_out.(subs{ii}).LTD1.close.trk_nsig;
        place.LTD1_trk(ii) = all_out.(subs{ii}).LTD1.place.trk_nsig;
        
        all.LTD1_post(ii) = all_out.(subs{ii}).LTD1.all.post_nsig;
        open.LTD1_post(ii) = all_out.(subs{ii}).LTD1.open.post_nsig;
        closed.LTD1_post(ii) = all_out.(subs{ii}).LTD1.close.post_nsig;
        place.LTD1_post(ii) = all_out.(subs{ii}).LTD1.place.post_nsig;
    end
    
    if isfield(all_out.(subs{ii}), 'LTD5')
        all.LTD5_pre(ii) = all_out.(subs{ii}).LTD5.all.pre_nsig;
        open.LTD5_pre(ii) = all_out.(subs{ii}).LTD5.open.pre_nsig;
        closed.LTD5_pre(ii) = all_out.(subs{ii}).LTD5.close.pre_nsig;
        place.LTD5_pre(ii) = all_out.(subs{ii}).LTD5.place.pre_nsig;
        
        all.LTD5_trk(ii) = all_out.(subs{ii}).LTD5.all.trk_nsig;
        open.LTD5_trk(ii) = all_out.(subs{ii}).LTD5.open.trk_nsig;
        closed.LTD5_trk(ii) = all_out.(subs{ii}).LTD5.close.trk_nsig;
        place.LTD5_trk(ii) = all_out.(subs{ii}).LTD5.place.trk_nsig;
        
        all.LTD5_post(ii) = all_out.(subs{ii}).LTD5.all.post_nsig;
        open.LTD5_post(ii) = all_out.(subs{ii}).LTD5.open.post_nsig;
        closed.LTD5_post(ii) = all_out.(subs{ii}).LTD5.close.post_nsig;
        place.LTD5_post(ii) = all_out.(subs{ii}).LTD5.place.post_nsig;
    end
    
    
    if isfield(all_out.(subs{ii}), 'HATD1')
        all.HATD1_pre(ii) = all_out.(subs{ii}).HATD1.all.pre_nsig;
        open.HATD1_pre(ii) = all_out.(subs{ii}).HATD1.open.pre_nsig;
        closed.HATD1_pre(ii) = all_out.(subs{ii}).HATD1.close.pre_nsig;
        place.HATD1_pre(ii) = all_out.(subs{ii}).HATD1.place.pre_nsig;
        
        all.HATD1_trk(ii) = all_out.(subs{ii}).HATD1.all.trk_nsig;
        open.HATD1_trk(ii) = all_out.(subs{ii}).HATD1.open.trk_nsig;
        closed.HATD1_trk(ii) = all_out.(subs{ii}).HATD1.close.trk_nsig;
        place.HATD1_trk(ii) = all_out.(subs{ii}).HATD1.place.trk_nsig;
        
        all.HATD1_post(ii) = all_out.(subs{ii}).HATD1.all.post_nsig;
        open.HATD1_post(ii) = all_out.(subs{ii}).HATD1.open.post_nsig;
        closed.HATD1_post(ii) = all_out.(subs{ii}).HATD1.close.post_nsig;
        place.HATD1_post(ii) = all_out.(subs{ii}).HATD1.place.post_nsig;
    end
    
    if isfield(all_out.(subs{ii}), 'HATD5')
        all.HATD5_pre(ii) = all_out.(subs{ii}).HATD5.all.pre_nsig;
        open.HATD5_pre(ii) = all_out.(subs{ii}).HATD5.open.pre_nsig;
        closed.HATD5_pre(ii) = all_out.(subs{ii}).HATD5.close.pre_nsig;
        place.HATD5_pre(ii) = all_out.(subs{ii}).HATD5.place.pre_nsig;
        
        all.HATD5_trk(ii) = all_out.(subs{ii}).HATD5.all.trk_nsig;
        open.HATD5_trk(ii) = all_out.(subs{ii}).HATD5.open.trk_nsig;
        closed.HATD5_trk(ii) = all_out.(subs{ii}).HATD5.close.trk_nsig;
        place.HATD5_trk(ii) = all_out.(subs{ii}).HATD5.place.trk_nsig;
        
        all.HATD5_post(ii) = all_out.(subs{ii}).HATD5.all.post_nsig;
        open.HATD5_post(ii) = all_out.(subs{ii}).HATD5.open.post_nsig;
        closed.HATD5_post(ii) = all_out.(subs{ii}).HATD5.close.post_nsig;
        place.HATD5_post(ii) = all_out.(subs{ii}).HATD5.place.post_nsig;
    end
    
    if isfield(all_out.(subs{ii}), 'HATDSwitch')
        all.HATDS_pre(ii) = all_out.(subs{ii}).HATDSwitch.all.pre_nsig;
        open.HATDS_pre(ii) = all_out.(subs{ii}).HATDSwitch.open.pre_nsig;
        closed.HATDS_pre(ii) = all_out.(subs{ii}).HATDSwitch.close.pre_nsig;
        place.HATDS_pre(ii) = all_out.(subs{ii}).HATDSwitch.place.pre_nsig;
        
        all.HATDS_trk(ii) = all_out.(subs{ii}).HATDSwitch.all.trk_nsig;
        open.HATDS_trk(ii) = all_out.(subs{ii}).HATDSwitch.open.trk_nsig;
        closed.HATDS_trk(ii) = all_out.(subs{ii}).HATDSwitch.close.trk_nsig;
        place.HATDS_trk(ii) = all_out.(subs{ii}).HATDSwitch.place.trk_nsig;
        
        all.HATDS_post(ii) = all_out.(subs{ii}).HATDSwitch.all.post_nsig;
        open.HATDS_post(ii) = all_out.(subs{ii}).HATDSwitch.open.post_nsig;
        closed.HATDS_post(ii) = all_out.(subs{ii}).HATDSwitch.close.post_nsig;
        place.HATDS_post(ii) = all_out.(subs{ii}).HATDSwitch.place.post_nsig;
    end
end

%% plot bars 
figure(301)
clf
plot_type = fieldnames(all);
m = 4; 
n = 5; 
c_ord = linspecer(40); 

% all cells
subplot(m, n, 1)
vals = [all.LTD1_pre, all.LTD1_trk, all.LTD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(1,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD1 all cells')
xlim([.25 3.75])
ylabel('% sig corr pairs');

subplot(m, n, 2)
vals = [all.LTD5_pre, all.LTD5_trk, all.LTD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(2,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD5 all cells')
    xlim([.25 3.75])

subplot(m, n, 3)
vals = [all.HATD1_pre, all.HATD1_trk, all.HATD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(3,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))     
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD1 all cells')
xlim([.25 3.75])

subplot(m, n, 4)
vals = [all.HATD5_pre, all.HATD5_trk, all.HATD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(4,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD5 all cells');
xlim([.25 3.75])

subplot(m, n, 5)
vals = [all.HATDS_pre, all.HATDS_trk, all.HATDS_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(5,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATDS all cells');
xlim([.25 3.75])


% all cells
subplot(m, n, 6)
vals = [open.LTD1_pre, open.LTD1_trk, open.LTD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(16,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD1 open cells')
xlim([.25 3.75])
ylabel('% sig corr pairs');

subplot(m, n, 7)
vals = [open.LTD5_pre, open.LTD5_trk, open.LTD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(17,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD5 open cells')
    xlim([.25 3.75])

subplot(m, n, 8)
vals = [open.HATD1_pre, open.HATD1_trk, open.HATD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(18,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))     
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD1 open cells')
xlim([.25 3.75])

subplot(m, n, 9)
vals = [open.HATD5_pre, open.HATD5_trk, open.HATD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(19,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD5 open cells');
xlim([.25 3.75])

subplot(m, n, 10)
vals = [open.HATDS_pre, open.HATDS_trk, open.HATDS_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(20,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATDS open cells');
xlim([.25 3.75])



% closed cells
subplot(m, n, 11)
vals = [closed.LTD1_pre, closed.LTD1_trk, closed.LTD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(8,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD1 closed cells')
xlim([.25 3.75])
ylabel('% sig corr pairs');

subplot(m, n, 12)
vals = [closed.LTD5_pre, closed.LTD5_trk, closed.LTD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(9,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD5 closed cells')
    xlim([.25 3.75])

subplot(m, n, 13)
vals = [closed.HATD1_pre, closed.HATD1_trk, closed.HATD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(10,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))     
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD1 closed cells')
xlim([.25 3.75])

subplot(m, n, 14)
vals = [closed.HATD5_pre, closed.HATD5_trk, closed.HATD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(11,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD5 closed cells');
xlim([.25 3.75])

subplot(m, n, 15)
vals = [closed.HATDS_pre, closed.HATDS_trk, closed.HATDS_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(12,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATDS closed cells');
xlim([.25 3.75])


% place cells
subplot(m, n, 16)
vals = [place.LTD1_pre, place.LTD1_trk, place.LTD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(36,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD1 place cells')
xlim([.25 3.75])
ylabel('% sig corr pairs');

subplot(m, n, 17)
vals = [place.LTD5_pre, place.LTD5_trk, place.LTD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(37,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD5 closed cells')
    xlim([.25 3.75])

subplot(m, n, 18)
vals = [place.HATD1_pre, place.HATD1_trk, place.HATD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(38,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))     
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD1 place cells')
xlim([.25 3.75])

subplot(m, n, 19)
vals = [place.HATD5_pre, place.HATD5_trk, place.HATD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(39,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD5 place cells');
xlim([.25 3.75])

subplot(m, n, 20)
vals = [place.HATDS_pre, place.HATDS_trk, place.HATDS_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(40,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATDS place cells');
xlim([.25 3.75])

saveas(gcf, [inter_dir filesep 'Summary_Sig_corr.fig'])
saveas(gcf, [inter_dir filesep 'Summary_Sig_corr.png'])


%% Get the overlap between sig correlated cells across phases

subs = fieldnames(all_out);

all.LTD1_pre = NaN(size(subs)); all.LTD1_trk = NaN(size(subs)); all.LTD1_post = NaN(size(subs));
all.LTD5_pre = NaN(size(subs)); all.LTD5_trk = NaN(size(subs)); all.LTD5_post = NaN(size(subs));
all.HATD1_pre = NaN(size(subs)); all.HATD1_trk = NaN(size(subs)); all.HATD1_post = NaN(size(subs));
all.HATD5_pre = NaN(size(subs)); all.HATD5_trk = NaN(size(subs)); all.HATD5_post = NaN(size(subs));
all.HATDS_pre = NaN(size(subs)); all.HATDS_trk = NaN(size(subs)); all.HATDS_post = NaN(size(subs));
open = all; closed = all; place = all;

for ii = 1:length(subs)
    if isfield(all_out.(subs{ii}), 'LTD1')
        all.LTD1_pre_trk(ii) = nnz(all_out.(subs{ii}).LTD1.all.pre_sig_ind & all_out.(subs{ii}).LTD1.all.trk_sig_ind)/numel(all_out.(subs{ii}).LTD1.all.pre_sig_ind); 
        open.LTD1_pre(ii) = all_out.(subs{ii}).LTD1.open.pre_nsig;
        closed.LTD1_pre(ii) = all_out.(subs{ii}).LTD1.close.pre_nsig;
        place.LTD1_pre(ii) = all_out.(subs{ii}).LTD1.place.pre_nsig;
        
        all.LTD1_trk(ii) = all_out.(subs{ii}).LTD1.all.trk_nsig;
        open.LTD1_trk(ii) = all_out.(subs{ii}).LTD1.open.trk_nsig;
        closed.LTD1_trk(ii) = all_out.(subs{ii}).LTD1.close.trk_nsig;
        place.LTD1_trk(ii) = all_out.(subs{ii}).LTD1.place.trk_nsig;
        
        all.LTD1_post(ii) = all_out.(subs{ii}).LTD1.all.post_nsig;
        open.LTD1_post(ii) = all_out.(subs{ii}).LTD1.open.post_nsig;
        closed.LTD1_post(ii) = all_out.(subs{ii}).LTD1.close.post_nsig;
        place.LTD1_post(ii) = all_out.(subs{ii}).LTD1.place.post_nsig;
    end
    
    if isfield(all_out.(subs{ii}), 'LTD5')
        all.LTD5_pre(ii) = all_out.(subs{ii}).LTD5.all.pre_nsig;
        open.LTD5_pre(ii) = all_out.(subs{ii}).LTD5.open.pre_nsig;
        closed.LTD5_pre(ii) = all_out.(subs{ii}).LTD5.close.pre_nsig;
        place.LTD5_pre(ii) = all_out.(subs{ii}).LTD5.place.pre_nsig;
        
        all.LTD5_trk(ii) = all_out.(subs{ii}).LTD5.all.trk_nsig;
        open.LTD5_trk(ii) = all_out.(subs{ii}).LTD5.open.trk_nsig;
        closed.LTD5_trk(ii) = all_out.(subs{ii}).LTD5.close.trk_nsig;
        place.LTD5_trk(ii) = all_out.(subs{ii}).LTD5.place.trk_nsig;
        
        all.LTD5_post(ii) = all_out.(subs{ii}).LTD5.all.post_nsig;
        open.LTD5_post(ii) = all_out.(subs{ii}).LTD5.open.post_nsig;
        closed.LTD5_post(ii) = all_out.(subs{ii}).LTD5.close.post_nsig;
        place.LTD5_post(ii) = all_out.(subs{ii}).LTD5.place.post_nsig;
    end
    
    
    if isfield(all_out.(subs{ii}), 'HATD1')
        all.HATD1_pre(ii) = all_out.(subs{ii}).HATD1.all.pre_nsig;
        open.HATD1_pre(ii) = all_out.(subs{ii}).HATD1.open.pre_nsig;
        closed.HATD1_pre(ii) = all_out.(subs{ii}).HATD1.close.pre_nsig;
        place.HATD1_pre(ii) = all_out.(subs{ii}).HATD1.place.pre_nsig;
        
        all.HATD1_trk(ii) = all_out.(subs{ii}).HATD1.all.trk_nsig;
        open.HATD1_trk(ii) = all_out.(subs{ii}).HATD1.open.trk_nsig;
        closed.HATD1_trk(ii) = all_out.(subs{ii}).HATD1.close.trk_nsig;
        place.HATD1_trk(ii) = all_out.(subs{ii}).HATD1.place.trk_nsig;
        
        all.HATD1_post(ii) = all_out.(subs{ii}).HATD1.all.post_nsig;
        open.HATD1_post(ii) = all_out.(subs{ii}).HATD1.open.post_nsig;
        closed.HATD1_post(ii) = all_out.(subs{ii}).HATD1.close.post_nsig;
        place.HATD1_post(ii) = all_out.(subs{ii}).HATD1.place.post_nsig;
    end
    
    if isfield(all_out.(subs{ii}), 'HATD5')
        all.HATD5_pre(ii) = all_out.(subs{ii}).HATD5.all.pre_nsig;
        open.HATD5_pre(ii) = all_out.(subs{ii}).HATD5.open.pre_nsig;
        closed.HATD5_pre(ii) = all_out.(subs{ii}).HATD5.close.pre_nsig;
        place.HATD5_pre(ii) = all_out.(subs{ii}).HATD5.place.pre_nsig;
        
        all.HATD5_trk(ii) = all_out.(subs{ii}).HATD5.all.trk_nsig;
        open.HATD5_trk(ii) = all_out.(subs{ii}).HATD5.open.trk_nsig;
        closed.HATD5_trk(ii) = all_out.(subs{ii}).HATD5.close.trk_nsig;
        place.HATD5_trk(ii) = all_out.(subs{ii}).HATD5.place.trk_nsig;
        
        all.HATD5_post(ii) = all_out.(subs{ii}).HATD5.all.post_nsig;
        open.HATD5_post(ii) = all_out.(subs{ii}).HATD5.open.post_nsig;
        closed.HATD5_post(ii) = all_out.(subs{ii}).HATD5.close.post_nsig;
        place.HATD5_post(ii) = all_out.(subs{ii}).HATD5.place.post_nsig;
    end
    
    if isfield(all_out.(subs{ii}), 'HATDSwitch')
        all.HATDS_pre(ii) = all_out.(subs{ii}).HATDSwitch.all.pre_nsig;
        open.HATDS_pre(ii) = all_out.(subs{ii}).HATDSwitch.open.pre_nsig;
        closed.HATDS_pre(ii) = all_out.(subs{ii}).HATDSwitch.close.pre_nsig;
        place.HATDS_pre(ii) = all_out.(subs{ii}).HATDSwitch.place.pre_nsig;
        
        all.HATDS_trk(ii) = all_out.(subs{ii}).HATDSwitch.all.trk_nsig;
        open.HATDS_trk(ii) = all_out.(subs{ii}).HATDSwitch.open.trk_nsig;
        closed.HATDS_trk(ii) = all_out.(subs{ii}).HATDSwitch.close.trk_nsig;
        place.HATDS_trk(ii) = all_out.(subs{ii}).HATDSwitch.place.trk_nsig;
        
        all.HATDS_post(ii) = all_out.(subs{ii}).HATDSwitch.all.post_nsig;
        open.HATDS_post(ii) = all_out.(subs{ii}).HATDSwitch.open.post_nsig;
        closed.HATDS_post(ii) = all_out.(subs{ii}).HATDSwitch.close.post_nsig;
        place.HATDS_post(ii) = all_out.(subs{ii}).HATDSwitch.place.post_nsig;
    end
end

%% plot bars 
figure(301)
clf
plot_type = fieldnames(all);
m = 4; 
n = 5; 
c_ord = linspecer(40); 

% all cells
subplot(m, n, 1)
vals = [all.LTD1_pre, all.LTD1_trk, all.LTD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(1,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD1 all cells')
xlim([.25 3.75])
ylabel('% sig corr pairs');

subplot(m, n, 2)
vals = [all.LTD5_pre, all.LTD5_trk, all.LTD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(2,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD5 all cells')
    xlim([.25 3.75])

subplot(m, n, 3)
vals = [all.HATD1_pre, all.HATD1_trk, all.HATD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(3,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))     
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD1 all cells')
xlim([.25 3.75])

subplot(m, n, 4)
vals = [all.HATD5_pre, all.HATD5_trk, all.HATD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(4,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD5 all cells');
xlim([.25 3.75])

subplot(m, n, 5)
vals = [all.HATDS_pre, all.HATDS_trk, all.HATDS_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(5,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATDS all cells');
xlim([.25 3.75])


% all cells
subplot(m, n, 6)
vals = [open.LTD1_pre, open.LTD1_trk, open.LTD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(16,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD1 open cells')
xlim([.25 3.75])
ylabel('% sig corr pairs');

subplot(m, n, 7)
vals = [open.LTD5_pre, open.LTD5_trk, open.LTD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(17,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD5 open cells')
    xlim([.25 3.75])

subplot(m, n, 8)
vals = [open.HATD1_pre, open.HATD1_trk, open.HATD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(18,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))     
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD1 open cells')
xlim([.25 3.75])

subplot(m, n, 9)
vals = [open.HATD5_pre, open.HATD5_trk, open.HATD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(19,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD5 open cells');
xlim([.25 3.75])

subplot(m, n, 10)
vals = [open.HATDS_pre, open.HATDS_trk, open.HATDS_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(20,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATDS open cells');
xlim([.25 3.75])



% closed cells
subplot(m, n, 11)
vals = [closed.LTD1_pre, closed.LTD1_trk, closed.LTD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(8,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD1 closed cells')
xlim([.25 3.75])
ylabel('% sig corr pairs');

subplot(m, n, 12)
vals = [closed.LTD5_pre, closed.LTD5_trk, closed.LTD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(9,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD5 closed cells')
    xlim([.25 3.75])

subplot(m, n, 13)
vals = [closed.HATD1_pre, closed.HATD1_trk, closed.HATD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(10,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))     
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD1 closed cells')
xlim([.25 3.75])

subplot(m, n, 14)
vals = [closed.HATD5_pre, closed.HATD5_trk, closed.HATD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(11,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD5 closed cells');
xlim([.25 3.75])

subplot(m, n, 15)
vals = [closed.HATDS_pre, closed.HATDS_trk, closed.HATDS_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(12,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATDS closed cells');
xlim([.25 3.75])


% place cells
subplot(m, n, 16)
vals = [place.LTD1_pre, place.LTD1_trk, place.LTD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(36,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD1 place cells')
xlim([.25 3.75])
ylabel('% sig corr pairs');

subplot(m, n, 17)
vals = [place.LTD5_pre, place.LTD5_trk, place.LTD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(37,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('LTD5 closed cells')
    xlim([.25 3.75])

subplot(m, n, 18)
vals = [place.HATD1_pre, place.HATD1_trk, place.HATD1_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(38,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))     
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD1 place cells')
xlim([.25 3.75])

subplot(m, n, 19)
vals = [place.HATD5_pre, place.HATD5_trk, place.HATD5_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(39,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))       
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATD5 place cells');
xlim([.25 3.75])

subplot(m, n, 20)
vals = [place.HATDS_pre, place.HATDS_trk, place.HATDS_post];
hold on
bar(1:3, nanmean(vals), 'facecolor', c_ord(40,:))
er = errorbar(1:3, nanmean(vals), nanstd(vals), 'color', 'k'); %./sqrt(sum(~isnan(nansum(vals,2))))    
set(gca,'xtick', 1:3,  'xticklabel', {'pre', 'trk', 'post'}); 
er.LineStyle =  'none';  
title('HATDS place cells');
xlim([.25 3.75])

saveas(gcf, [inter_dir filesep 'Summary_Sig_corr.fig'])
saveas(gcf, [inter_dir filesep 'Summary_Sig_corr.png'])