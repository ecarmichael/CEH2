%% sandbox_PCA_ICA_crawl

cd('C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter')

f_list = dir('pv*');

out_fig = [];
for ii =1:length(f_list)
    
    [out_fig.REM_z{ii}, out_fig.Az{ii}] = sandbox_PCA_ICA(f_list(ii).name);
    
    clearvars -except f_list out ii
    close all
end

%%  same thing but stats only


cd('C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter')
fig_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\Comp_Can_inter\Assembly\checks'; 
% cd('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter')

f_list = dir('*data*');



out = [];
session = []; novel_idx = []; anx_idx = []; HS_idx = [];
for ii = 1:length(f_list)
    session{ii} = f_list(ii).name;
        [A_out{ii}] = MS_PCA_ICA_no_fig(f_list(ii).name, fig_dir);
        close all
    
    if ~isempty(strfind(f_list(ii).name, 'HATDS'))
        HS_idx(ii) = 1;
    else
        HS_idx(ii) = 0;
    end
    
    if ~isempty(strfind(f_list(ii).name, 'D1')) %|| ~isempty(strfind(f_list(ii).name, 'HATDS'))
        novel_idx(ii) = 1;
    end
    
    if ~isempty(strfind(f_list(ii).name, 'D5'))
        novel_idx(ii) = 0;
    end
    
    if ~isempty(strfind(f_list(ii).name, 'HAT'))
        anx_idx(ii) = 1;
    else
        anx_idx(ii) = 0;
    end
    
    
    %    clearvars -except f_list out ii
    %    close all
end


%% same thing for J20 mice
data_dir = ('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV'); 
fig_dir = 'C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Assembly_EV\checks'; 
% cd('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter')
cd(data_dir); 
j20_list = dir('*data*');



out = [];
session = []; novel_idx = []; D3_idx = []; HS_idx = [];
for ii = 1:length(j20_list)
    session{ii} = j20_list(ii).name;
    cd(data_dir)
        [A_out{ii}] = MS_PCA_ICA_no_fig(j20_list(ii).name, fig_dir);
        close all
    

    
    if ~isempty(strfind(j20_list(ii).name, 'D1')) %|| ~isempty(strfind(f_list(ii).name, 'HATDS'))
        novel_idx(ii) = 1;
        D3_idx(ii) = 0;
    end
    
    if ~isempty(strfind(j20_list(ii).name, 'D3'))
        D3_idx(ii) = 1;
        novel_idx(ii) = 0;
    end

    
    if ~isempty(strfind(j20_list(ii).name, 'D5'))
        novel_idx(ii) = 0;
        D3_idx(ii) = 0;
    end

    
    
    %    clearvars -except f_list out ii
    %    close all
end


%% collect output
c_ord = MS_linspecer(5); 

HS_idx = logical(HS_idx(1:length(A_out))); 
n_idx = logical(novel_idx(1:length(A_out))); 
a_idx = logical(anx_idx(1:length(A_out))); 

n_idx = n_idx & ~HS_idx; 
a_idx = a_idx & ~HS_idx; 

lt1_idx = n_idx & ~a_idx & ~HS_idx;  
lt5_idx = ~n_idx & ~a_idx & ~HS_idx;  

H1_idx = n_idx & a_idx & ~HS_idx;  
H5_idx = ~n_idx & a_idx & ~HS_idx;  



% get the length of the assemblies in each recording. 
for ii =length(A_out):-1:1
    A_l(ii) = length(A_out{ii}.ReAct_rate_pre);
    W_l(ii) = size(A_out{ii}.Pos_templates,2);
end

Pre_A_rate = NaN(length(A_out), max(A_l)); 
Post_A_rate = NaN(length(A_out), max(A_l)); 
Wake_A_rate = NaN(length(A_out), max(W_l)); 

Pre_A_sig = []; Post_A_sig = [];
C_pre_idx = []; C_post_idx = [];
O_pre_idx = []; O_post_idx = [];
M_pre_idx = []; M_pre_idx = [];

for ii = length(A_out):-1:1
    pre_sig_idx = logical(A_out{ii}.Ass_p_val_pre < 0.05); 
    post_sig_idx = logical(A_out{ii}.Ass_p_val_post < 0.05); 

    Pre_A_sig(ii) = sum(pre_sig_idx);
    Post_A_sig(ii) = sum(post_sig_idx);
    Wake_A_sig(ii) = size(A_out{ii}.Pos_templates,2); 
    
    ReAct_sig(ii) = mean(A_out{ii}.ReAct_S(pre_sig_idx & post_sig_idx)); 
    
    Pre_A_rate(ii,1:sum(pre_sig_idx)) = A_out{ii}.ReAct_rate_pre(pre_sig_idx); 
    Post_A_rate(ii,1:sum(post_sig_idx)) = A_out{ii}.ReAct_rate_post(post_sig_idx);
    
    Wake_A_rate(ii,1:sum(pre_sig_idx)) = A_out{ii}.ReAct_rate_pre(pre_sig_idx); 

% %     Wake_A_rate(ii  
% for jj = size(A_out{ii}.Pos_projections,1):-1:1
%     Wake_A_Rate(ii,jj) = sum(A_out{ii}.Pos_projections(jj,:) > 8)/length(A_out{ii}.Pos_projections)/30/60; 
% end


    
    Pre_A_rate(ii,1:sum(pre_sig_idx)) = A_out{ii}.ReAct_rate_pre(pre_sig_idx); 
    Post_A_rate(ii,1:sum(post_sig_idx)) = A_out{ii}.ReAct_rate_post(post_sig_idx);
    
    C_pre_idx(ii) = sum(ismember(find(pre_sig_idx),A_out{ii}.close_idx))/length(pre_sig_idx); 
    O_pre_idx(ii) = sum(ismember(find(pre_sig_idx),A_out{ii}.open_idx))/length(pre_sig_idx); 
    M_pre_idx(ii) = sum(ismember(find(pre_sig_idx),A_out{ii}.mid_idx))/length(pre_sig_idx); 
    
    C_post_idx(ii) = sum(ismember(find(post_sig_idx),A_out{ii}.close_idx))/length(post_sig_idx); 
    O_post_idx(ii) = sum(ismember(find(post_sig_idx),A_out{ii}.open_idx))/length(post_sig_idx); 
    M_post_idx(ii) = sum(ismember(find(post_sig_idx),A_out{ii}.mid_idx))/length(post_sig_idx); 
end

n = 3; m = 5; 
figure(1010)
subplot(n,m,1)
cla
bar([mean(Pre_A_sig(n_idx)), mean(Post_A_sig(n_idx))], 'FaceColor', c_ord(1,:), 'EdgeColor', c_ord(1,:));
hold on
eb = errorbar([mean(Pre_A_sig(n_idx)), mean(Post_A_sig(n_idx))], [MS_SEM(Pre_A_sig(n_idx)) ,MS_SEM(Post_A_sig(n_idx))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Novel')
ylabel('# assemblies')


subplot(n,m,2)
cla
bar([mean(Pre_A_sig(~n_idx)), mean(Post_A_sig(~n_idx))], 'FaceColor', c_ord(4,:), 'EdgeColor', c_ord(4,:));
hold on
eb = errorbar([mean(Pre_A_sig(~n_idx)), mean(Post_A_sig(~n_idx))], [MS_SEM(Pre_A_sig(~n_idx)) ,MS_SEM(Post_A_sig(~n_idx))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Familiar')
% ylabel('# assemblies')

subplot(n,m,3)
cla
bar([mean(Pre_A_sig(~a_idx)), mean(Post_A_sig(~a_idx))], 'FaceColor', c_ord(3,:), 'EdgeColor', c_ord(3,:));
hold on
eb = errorbar([mean(Pre_A_sig(~a_idx)), mean(Post_A_sig(~a_idx))], [MS_SEM(Pre_A_sig(~a_idx)) ,MS_SEM(Post_A_sig(~a_idx))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Linear Track')

subplot(n,m,4)
cla
bar([mean(Pre_A_sig(a_idx)), mean(Post_A_sig(a_idx))], 'FaceColor', c_ord(2,:), 'EdgeColor', c_ord(2,:));
hold on
eb = errorbar([mean(Pre_A_sig(a_idx)), mean(Post_A_sig(a_idx))], [MS_SEM(Pre_A_sig(a_idx)) ,MS_SEM(Post_A_sig(a_idx))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Half-Anxiety')

subplot(n,m,5)
cla
bar([mean(Pre_A_sig(HS_idx)), mean(Post_A_sig(HS_idx))], 'FaceColor', c_ord(5,:), 'EdgeColor', c_ord(5,:));
hold on
eb = errorbar([mean(Pre_A_sig(HS_idx)), mean(Post_A_sig(HS_idx))], [MS_SEM(Pre_A_sig(HS_idx)) ,MS_SEM(Post_A_sig(HS_idx))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Half-Anxiety Switch')


% reactivation rate plots
subplot(n,m,6)
cla
bar([nanmean(Pre_A_rate(n_idx,:), 'all'), nanmean(Post_A_rate(n_idx,:), 'all')], 'FaceColor', c_ord(1,:), 'EdgeColor', c_ord(1,:));
hold on
eb = errorbar([nanmean(Pre_A_rate(n_idx,:), 'all'), nanmean(Post_A_rate(n_idx,:), 'all')], [MS_SEM(Pre_A_rate(n_idx)) ,MS_SEM(Post_A_rate(n_idx))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Novel')
ylabel('ReAct rate /min)')


subplot(n,m,7)
cla
bar([nanmean(Pre_A_rate(~n_idx,:), 'all'), nanmean(Post_A_rate(~n_idx,:), 'all')], 'FaceColor', c_ord(4,:), 'EdgeColor', c_ord(4,:));
hold on
eb = errorbar([nanmean(Pre_A_rate(~n_idx,:), 'all'), nanmean(Post_A_rate(~n_idx,:), 'all')], [MS_SEM(Pre_A_rate(~n_idx)) ,MS_SEM(Post_A_rate(~n_idx))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Familiar')


subplot(n,m,8)
cla
bar([nanmean(Pre_A_rate(a_idx,:), 'all'), nanmean(Post_A_rate(a_idx,:), 'all')], 'FaceColor', c_ord(3,:), 'EdgeColor', c_ord(3,:));
hold on
eb = errorbar([nanmean(Pre_A_rate(a_idx,:), 'all'), nanmean(Post_A_rate(a_idx,:), 'all')], [MS_SEM(Pre_A_rate(a_idx)) ,MS_SEM(Post_A_rate(a_idx))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Linear Track')


subplot(n,m,9)
cla
bar([nanmean(Pre_A_rate(~a_idx,:), 'all'), nanmean(Post_A_rate(~a_idx,:), 'all')], 'FaceColor', c_ord(2,:), 'EdgeColor', c_ord(2,:));
hold on
eb = errorbar([nanmean(Pre_A_rate(~a_idx,:), 'all'), nanmean(Post_A_rate(~a_idx,:), 'all')], [MS_SEM(Pre_A_rate(~a_idx)) ,MS_SEM(Post_A_rate(~a_idx))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Half-Anxiety')

subplot(n,m,10)
cla
bar([nanmean(Pre_A_rate(HS_idx)), nanmean(Post_A_rate(HS_idx))], 'FaceColor', c_ord(5,:), 'EdgeColor', c_ord(5,:));
hold on
eb = errorbar([nanmean(Pre_A_rate(HS_idx)), nanmean(Post_A_rate(HS_idx))], [MS_SEM(Pre_A_rate(HS_idx)) ,MS_SEM(Post_A_rate(HS_idx))]);
eb.LineStyle = 'none';
eb.Color = 'k';
set(gca,'xtick', 1:2, 'XTickLabel', {'Pre', 'Post'}, 'XTickLabelRotation', 45)
title('Half-Anxiety Switch')


%% per subject


figure(1011)
subplot(2,2,1)
cla
b = bar([nanmean(Pre_A_sig(lt1_idx)),nanmean(Wake_A_sig(lt1_idx)), nanmean(Post_A_sig(lt1_idx));...
    nanmean(Pre_A_sig(lt5_idx)),nanmean(Wake_A_sig(lt5_idx)), nanmean(Post_A_sig(lt5_idx));...
    nanmean(Pre_A_sig(H1_idx)),nanmean(Wake_A_sig(H1_idx)), nanmean(Post_A_sig(H1_idx));...
    nanmean(Pre_A_sig(H5_idx)),nanmean(Wake_A_sig(H5_idx)), nanmean(Post_A_sig(H5_idx));...
    nanmean(Pre_A_sig(HS_idx)),nanmean(Wake_A_sig(HS_idx)), nanmean(Post_A_sig(HS_idx));...
    ]);
% hold on
[~, h, eb] = errorbar_groups([nanmean(Pre_A_sig(lt1_idx)),nanmean(Wake_A_sig(lt1_idx)), nanmean(Post_A_sig(lt1_idx));...
    nanmean(Pre_A_sig(lt5_idx)),nanmean(Wake_A_sig(lt5_idx)), nanmean(Post_A_sig(lt5_idx));...
    nanmean(Pre_A_sig(H1_idx)),nanmean(Wake_A_sig(H1_idx)), nanmean(Post_A_sig(H1_idx));...
    nanmean(Pre_A_sig(H5_idx)),nanmean(Wake_A_sig(H5_idx)), nanmean(Post_A_sig(H5_idx));...
    nanmean(Pre_A_sig(HS_idx)),nanmean(Wake_A_sig(HS_idx)), nanmean(Post_A_sig(HS_idx))]',[...
    MS_SEM(Pre_A_sig(lt1_idx)),MS_SEM(Wake_A_sig(lt1_idx)), MS_SEM(Post_A_sig(lt1_idx));...
    MS_SEM(Pre_A_sig(lt5_idx)),MS_SEM(Wake_A_sig(lt5_idx)), MS_SEM(Post_A_sig(lt5_idx));...
    MS_SEM(Pre_A_sig(H1_idx)),MS_SEM(Wake_A_sig(H1_idx)), MS_SEM(Post_A_sig(H1_idx));...
    MS_SEM(Pre_A_sig(H5_idx)),MS_SEM(Wake_A_sig(H5_idx)), MS_SEM(Post_A_sig(H5_idx));...
    MS_SEM(Pre_A_sig(HS_idx)),MS_SEM(Wake_A_sig(HS_idx)), MS_SEM(Post_A_sig(HS_idx));...
    ]', 'FigID', gcf, 'AxID', gca);

ylabel('# sig assemblies')
set(gca, 'xticklabel', {'Novel', 'Fam.', 'Anxiety D1', 'Anxiety D5', 'Anxiety Switch'})
legend({'Pre', 'Wake', 'Post'})

%s Sig Ass Rate
subplot(2,2,1)
cla
b = bar([nanmean(Pre_A_sig(lt1_idx)),nanmean(Wake_A_sig(lt1_idx)), nanmean(Post_A_sig(lt1_idx));...
    nanmean(Pre_A_sig(lt5_idx)),nanmean(Wake_A_sig(lt5_idx)), nanmean(Post_A_sig(lt5_idx));...
    nanmean(Pre_A_sig(H1_idx)),nanmean(Wake_A_sig(H1_idx)), nanmean(Post_A_sig(H1_idx));...
    nanmean(Pre_A_sig(H5_idx)),nanmean(Wake_A_sig(H5_idx)), nanmean(Post_A_sig(H5_idx));...
    nanmean(Pre_A_sig(HS_idx)),nanmean(Wake_A_sig(HS_idx)), nanmean(Post_A_sig(HS_idx));...
    ]);
% hold on
[~, h, eb] = errorbar_groups([nanmean(Pre_A_sig(lt1_idx)),nanmean(Wake_A_sig(lt1_idx)), nanmean(Post_A_sig(lt1_idx));...
    nanmean(Pre_A_sig(lt5_idx)),nanmean(Wake_A_sig(lt5_idx)), nanmean(Post_A_sig(lt5_idx));...
    nanmean(Pre_A_sig(H1_idx)),nanmean(Wake_A_sig(H1_idx)), nanmean(Post_A_sig(H1_idx));...
    nanmean(Pre_A_sig(H5_idx)),nanmean(Wake_A_sig(H5_idx)), nanmean(Post_A_sig(H5_idx));...
    nanmean(Pre_A_sig(HS_idx)),nanmean(Wake_A_sig(HS_idx)), nanmean(Post_A_sig(HS_idx))]',[...
    MS_SEM(Pre_A_sig(lt1_idx)),MS_SEM(Wake_A_sig(lt1_idx)), MS_SEM(Post_A_sig(lt1_idx));...
    MS_SEM(Pre_A_sig(lt5_idx)),MS_SEM(Wake_A_sig(lt5_idx)), MS_SEM(Post_A_sig(lt5_idx));...
    MS_SEM(Pre_A_sig(H1_idx)),MS_SEM(Wake_A_sig(H1_idx)), MS_SEM(Post_A_sig(H1_idx));...
    MS_SEM(Pre_A_sig(H5_idx)),MS_SEM(Wake_A_sig(H5_idx)), MS_SEM(Post_A_sig(H5_idx));...
    MS_SEM(Pre_A_sig(HS_idx)),MS_SEM(Wake_A_sig(HS_idx)), MS_SEM(Post_A_sig(HS_idx));...
    ]', 'FigID', gcf, 'AxID', gca);

ylabel('# sig assemblies')
set(gca, 'xticklabel', {'Novel', 'Fam.', 'Anxiety D1', 'Anxiety D5', 'Anxiety Switch'})
legend({'Pre', 'Wake', 'Post'})

% bar(
%% collect the output values
warning off
var_name = {'session', 'novel', 'anxiety', 'nAssemblies_z', 'nPlace_Assemblies','nREM_Assemblies', 'nWake_open', 'nWake_close','Wake_OC_idx', 'nREM_open', 'nREM_close','nREM_mid', 'REM_OC_idx'};

for ii = length(session):-1:1
    
    
    

    
    nAssemblies_z = cell2mat(out.Az);
    
    if isstruct(out.REM_z{ii})  %&& isempty(strfind(session{ii}, '12')) && isempty(strfind(session{ii}, '11'))
        nPlace_Assemblies(ii) = size(out.REM_z{ii}.all,1);
        nWake_open(ii) = sum(out.REM_z{ii}.isopen);
        nWake_close(ii) = sum(~out.REM_z{ii}.isopen);
        
        
        for iS = length(out.REM_z{ii}.shuff_time_prog_rem_z):-1:1
            all_z_peaks = [];
            for jj = size(out.REM_z{ii}.shuff_time_prog_rem_z{iS},1):-1:1
                [~, idx] = findpeaks(out.REM_z{ii}.shuff_time_prog_rem_z{iS}(jj,:), 'MinPeakHeight', 2, 'MinPeakDistance', 10);
                all_z_peaks(jj) = length(idx);
                
            end
            REM_shuff(ii,iS) = sum(all_z_peaks >0);
        end
        
        
        % All wake assemblies detected and replayed during REM. 
        all_a_peaks = [];
        for jj = size(out.REM_z{ii}.all_time_prog,1):-1:1
            [~, idx] = findpeaks(out.REM_z{ii}.all_time_prog_z(jj,:), 'MinPeakHeight', 5, 'MinPeakDistance', 10);
            all_a_peaks(jj) = length(idx);
            
        end
        nREM_all(ii) = (sum(all_a_peaks > 0) / size(out.REM_z{ii}.all_time_prog,1))*100; 
        nREM_all_raw(ii) = sum(all_a_peaks > 0) ; 

        open_peaks = []; W_open_peaks = [];
        for jj = size(out.REM_z{ii}.open,1):-1:1
            [~, idx] = findpeaks(out.REM_z{ii}.open(jj,:), 'MinPeakHeight', 5, 'MinPeakDistance', 10);
            open_peaks(jj) = length(idx);
            
%             this_z = zscore(out.wake{ii}(jj,:));
            [~, idx] = findpeaks(out.wake{ii}(jj,:), 'MinPeakHeight', 5, 'MinPeakDistance', 10);
            W_open_peaks(jj) = length(idx);
            
        end
        nREM_open(ii) = sum(open_peaks >0);
        nWake_open(ii) = sum(W_open_peaks >0);
        
        
        %mid assemblies
        mid_peaks = []; W_mid_peaks = [];
        for jj = size(out.REM_z{ii}.mid,1):-1:1
            [~, idx] = findpeaks(out.REM_z{ii}.mid(jj,:), 'MinPeakHeight', 5, 'MinPeakDistance', 10);
            mid_peaks(jj) = length(idx);
            
%             this_z = zscore();
            [~, idx] = findpeaks(out.wake{ii}(jj,:), 'MinPeakHeight', 5, 'MinPeakDistance', 10);
            W_mid_peaks(jj) = length(idx);
            
        end
        nREM_mid(ii) = sum(mid_peaks >0);
        nWake_mid(ii) = sum(W_mid_peaks >0);
        
        
        
        % same thing for closed REM events
        closed_peaks = []; W_closed_peaks = [];
        for jj = size(out.REM_z{ii}.close,1):-1:1
            [~, idx] = findpeaks(out.REM_z{ii}.close(jj,:), 'MinPeakHeight', 5, 'MinPeakDistance', 10);
            closed_peaks(jj) = length(idx);
            
            
            this_z = (out.wake{ii}(jj,:)); %zscore
            [~, idx] = findpeaks(this_z, 'MinPeakHeight', 5, 'MinPeakDistance', 10);
            W_closed_peaks(jj) = length(idx);
            
        end
        nREM_close(ii) = sum(closed_peaks > 0);
        nWake_close(ii) = sum(W_closed_peaks >0);
        
        Wake_OC_idx(ii) = (nWake_open(ii) - nWake_close(ii)) / (nWake_open(ii) + nWake_close(ii));
        REM_OC_idx(ii) = (nREM_open(ii) - nREM_close(ii)) / (nREM_open(ii) + nREM_mid(ii) + nREM_close(ii));
        REM_mid_OC_idx(ii) = (nWake_open(ii) - nWake_close(ii)) / (nREM_open(ii) - nREM_close(ii));

        
        nREM_Assemblies(ii) = nREM_open(ii) + nREM_close(ii);
        
    else
        nPlace_Assemblies(ii) = NaN;
        nREM_Assemblies(ii) = NaN;
        nREM_all_raw(ii) = NaN;
         nREM_all(ii) =NaN; 
        nWake_open(ii) = NaN;
        nWake_close(ii) = NaN;
        nWake_mid(ii) = NaN;
        nREM_close(ii) = NaN;
        nREM_open(ii) = NaN;
        nREM_mid(ii) = NaN; 
        REM_OC_idx(ii) = NaN;
        Wake_OC_idx(ii) = NaN;
        
    end
end

non_react_idx = (nREM_open == 0) | (nREM_close==0);
REM_OC_idx_act = REM_OC_idx;
REM_OC_idx_act(non_react_idx) = NaN;

pREM_open = (nREM_open./(nREM_open + nREM_mid + nREM_close))*100;
pREM_mid = (nREM_mid./(nREM_open + nREM_mid + nREM_close))*100;
pREM_close = (nREM_close./(nREM_open + nREM_mid + nREM_close))*100;


Ass_tbl = table(session, novel_idx, anx_idx, nAssemblies_z, nPlace_Assemblies,nREM_Assemblies, nWake_open, nWake_close,Wake_OC_idx,  nREM_open, nREM_close,nREM_mid, REM_OC_idx,'VariableNames', var_name );

anx_idx = logical(anx_idx);
novel_idx = logical(novel_idx);

fprintf('Chance level of Wake assemblies replayed in REM : %0.2f assemblies per session \n', nanmean(REM_shuff, 'all'))

%%  plot some stuff


HS_idx = logical(HS_idx);



p_ord = parula(9);

figure(909)
clf
subplot(2,3,1)
means_n_pA = [nanmean(nPlace_Assemblies(~anx_idx)); nanmean(nPlace_Assemblies(anx_idx)); nanmean(nPlace_Assemblies(novel_idx)); nanmean(nPlace_Assemblies(~novel_idx)), ; nanmean(nPlace_Assemblies(HS_idx))];
sem_n_pA = [MS_SEM(nPlace_Assemblies(~anx_idx)); MS_SEM(nPlace_Assemblies(anx_idx)); MS_SEM(nPlace_Assemblies(novel_idx)); MS_SEM(nPlace_Assemblies(~novel_idx));  MS_SEM(nPlace_Assemblies(HS_idx))];
hold on
%
eb = errorbar(1:5,means_n_pA, sem_n_pA);
eb.LineStyle = 'none';
b= bar(1:5,means_n_pA);
% boxplot([nPlace_Assemblies nPlace_Assemblies], [anx_idx   novel_idx+2])
b.FaceColor = p_ord(3,:);
b.EdgeColor = p_ord(3,:);
set(gca,'xtick', 1:5, 'XTickLabel', {'LT', 'HAT', 'Novel', 'Familiar', 'HATS'}, 'XTickLabelRotation', 45)
ylabel({'number of wake' ; 'place assemblies'})




subplot(2,3,2)

boxplot([Wake_OC_idx Wake_OC_idx Wake_OC_idx], [anx_idx   novel_idx+2 HS_idx+4])
set(gca,'xtick', 1:5, 'XTickLabel', {'LT', 'HAT', 'Novel', 'Familiar', 'HATS'}, 'XTickLabelRotation', 45, 'ytick', -1:1, 'YTickLabel', {'closed', '0', 'open'})
ylabel({'Wake assembly bias'})
xlim([.5 5.5])


subplot(2,3,4)
means_n_RA = [nanmean(nREM_Assemblies(~anx_idx)); nanmean(nREM_Assemblies(anx_idx)); nanmean(nREM_Assemblies(novel_idx)); nanmean(nREM_Assemblies(~novel_idx)); nanmean(nREM_Assemblies(HS_idx))];
sem_n_RA = [MS_SEM(nREM_Assemblies(~anx_idx)); MS_SEM(nREM_Assemblies(anx_idx)); MS_SEM(nREM_Assemblies(novel_idx)); MS_SEM(nREM_Assemblies(~novel_idx)); MS_SEM(nREM_Assemblies(HS_idx))];
hold on
%
eb = errorbar(1:5,means_n_RA, sem_n_RA);
eb.LineStyle = 'none';
b = bar(1:5,means_n_RA);
% boxplot([nPlace_Assemblies nPlace_Assemblies], [anx_idx   novel_idx+2])
b.FaceColor = p_ord(5,:);
b.EdgeColor = p_ord(5,:);
set(gca,'xtick', 1:5, 'XTickLabel', {'LT', 'HAT', 'Novel', 'Familiar', 'HATS'}, 'XTickLabelRotation', 45)
ylabel({'number of REM' ; 'reactivations'})




subplot(2,3,6)

means_n_RA = [nanmean(pREM_mid(~anx_idx)); nanmean(pREM_mid(anx_idx)); nanmean(pREM_mid(novel_idx)); nanmean(pREM_mid(~novel_idx)); nanmean(pREM_mid(HS_idx))];
sem_n_RA = [MS_SEM(pREM_mid(~anx_idx)); MS_SEM(pREM_mid(anx_idx)); MS_SEM(pREM_mid(novel_idx)); MS_SEM(pREM_mid(~novel_idx)); MS_SEM(pREM_mid(HS_idx))];
hold on
%
eb = errorbar(1:5,means_n_RA, sem_n_RA);
eb.LineStyle = 'none';
b = bar(1:5,means_n_RA);
% boxplot([nPlace_Assemblies nPlace_Assemblies], [anx_idx   novel_idx+2])
b.FaceColor = p_ord(5,:);
b.EdgeColor = p_ord(5,:);
set(gca,'xtick', 1:5, 'XTickLabel', {'LT', 'HAT', 'Novel', 'Familiar', 'HATS'}, 'XTickLabelRotation', 45)
ylabel({'REM % mid track reactivations'})
xlim([.5 5.5])




% subplot(2,3,5)
% % R_bias = [(REM_OC_idx(~anx_idx)), (REM_OC_idx(anx_idx)), (REM_OC_idx(novel_idx)), (REM_OC_idx(~novel_idx))];
% 
% boxplot([REM_OC_idx REM_OC_idx REM_OC_idx],  [anx_idx   novel_idx+2 HS_idx+4])
% 
% set(gca,'xtick', 1:5, 'XTickLabel', {'LT', 'HAT', 'Novel', 'Familiar', 'HATS'}, 'XTickLabelRotation', 45, 'ytick', -1:1, 'YTickLabel', {'closed', '0', 'open'})
% ylabel({'REM assembly bias'})
% xlim([.5 5.5])



%% simple figure
figure(102);
 clf
subplot(1,2,1)
means_n_pA = [nanmean(nREM_all_raw(novel_idx)); nanmean(nREM_all_raw(~novel_idx))];
sem_n_pA = [ MS_SEM(nREM_all_raw(novel_idx)); MS_SEM(nREM_all_raw(~novel_idx))];
hold on
%

eb = errorbar(1:2,means_n_pA', sem_n_pA);
eb.LineStyle = 'none';
eb.Color = 'k';

b= bar(1:2,means_n_pA');

% boxplot([nPlace_Assemblies nPlace_Assemblies], [anx_idx   novel_idx+2])
b.FaceColor = p_ord(3,:);
b.EdgeColor = p_ord(3,:);

means_n_pA = [nanmean(nREM_all_raw(anx_idx)); nanmean(nREM_all_raw(~anx_idx))];
sem_n_pA = [ MS_SEM(nREM_all_raw(anx_idx)); MS_SEM(nREM_all_raw(~anx_idx))];
hold on
%
eb = errorbar(3:4,means_n_pA', sem_n_pA);
eb.LineStyle = 'none';
eb.Color = 'k';
b= bar(3:4,means_n_pA');

% boxplot([nPlace_Assemblies nPlace_Assemblies], [anx_idx   novel_idx+2])
b.FaceColor = p_ord(8,:);
b.EdgeColor = p_ord(8,:);

set(gca,'xtick', 1:4, 'XTickLabel', {'Novel', 'Familiar', 'Linear', 'Anxiety'}, 'XTickLabelRotation', 45)
ylabel({'number of wake assemblies' ; ' Reactivated in REM '})
% xlim([.5 2.5])
axis square
[h, p] = ttest2(nREM_Assemblies(~novel_idx), nREM_Assemblies(novel_idx))
hline(4,'r', 'chance')
set(gca,'children',flipud(get(gca,'children')))

set(gca ,'Layer', 'Top')

subplot(1,2,2)
means_n_RA = [nanmean(pREM_mid(~anx_idx)); nanmean(pREM_mid(anx_idx))];
sem_n_RA = [MS_SEM(pREM_mid(~anx_idx)); MS_SEM(pREM_mid(anx_idx))];
hold on
%
eb = errorbar(1:2,means_n_RA, sem_n_RA);
eb.LineStyle = 'none';
eb.Color = 'k';

b = bar(1:2,means_n_RA);
% boxplot([nPlace_Assemblies nPlace_Assemblies], [anx_idx   novel_idx+2])
b.FaceColor = p_ord(6,:);
b.EdgeColor = p_ord(6,:);
set(gca,'xtick', 1:4, 'XTickLabel', {'Linear', 'Anxiety'}, 'XTickLabelRotation', 45)
ylabel({'REM % mid track reactivations'})
xlim([.5 2.5])
axis square
set(gcf, 'position', [1600 50 600 500])

% [h, p] = ttest2(pREM_mid(~anx_idx), pREM_mid(anx_idx))
% [p,t,stats] = anova1(nREM_mid, anx_idx);

% f = 
set(gca ,'Layer', 'Top')

% exportgraphics(gcf, ['C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\CIHR_2023_September\Assembly' filesep 'Assembly_Session_summary.pdf'], 'ContentType', 'vector')