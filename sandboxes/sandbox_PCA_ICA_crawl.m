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
% cd('/home/williamslab/Williams Lab Dropbox/Eric Carmichael/Comp_Can_inter')

f_list = dir('pv*');

out = [];
session = []; novel_idx = []; anx_idx = []; HS_idx = [];
for ii = 1:length(f_list)
    session{ii} = f_list(ii).name;
        [out.REM_z{ii}, out.Az{ii}, out.wake{ii}, out.rem{ii}] = sandbox_PCA_ICA_no_fig(f_list(ii).name);
    
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
            [~, idx] = findpeaks(out.REM_z{ii}.all_time_prog_z(jj,:), 'MinPeakHeight', 2, 'MinPeakDistance', 10);
            all_a_peaks(jj) = length(idx);
            
        end
        nREM_all(ii) = (sum(all_a_peaks > 0) / size(out.REM_z{ii}.all_time_prog,1))*100; 
        nREM_all_raw(ii) = sum(all_a_peaks > 0) ; 

        open_peaks = []; W_open_peaks = [];
        for jj = size(out.REM_z{ii}.open,1):-1:1
            [~, idx] = findpeaks(out.REM_z{ii}.open(jj,:), 'MinPeakHeight', 2, 'MinPeakDistance', 10);
            open_peaks(jj) = length(idx);
            
            this_z = zscore(out.wake{ii}(jj,:));
            [~, idx] = findpeaks(this_z, 'MinPeakHeight', 2, 'MinPeakDistance', 10);
            W_open_peaks(jj) = length(idx);
            
        end
        nREM_open(ii) = sum(open_peaks >0);
        nWake_open(ii) = sum(W_open_peaks >0);
        
        
        %mid assemblies
        mid_peaks = []; W_mid_peaks = [];
        for jj = size(out.REM_z{ii}.mid,1):-1:1
            [~, idx] = findpeaks(out.REM_z{ii}.mid(jj,:), 'MinPeakHeight', 2, 'MinPeakDistance', 10);
            mid_peaks(jj) = length(idx);
            
            this_z = zscore(out.wake{ii}(jj,:));
            [~, idx] = findpeaks(this_z, 'MinPeakHeight', 2, 'MinPeakDistance', 10);
            W_mid_peaks(jj) = length(idx);
            
        end
        nREM_mid(ii) = sum(mid_peaks >0);
        nWake_mid(ii) = sum(W_mid_peaks >0);
        
        
        
        % same thing for closed REM events
        closed_peaks = []; W_closed_peaks = [];
        for jj = size(out.REM_z{ii}.close,1):-1:1
            [~, idx] = findpeaks(out.REM_z{ii}.close(jj,:), 'MinPeakHeight', 1.96, 'MinPeakDistance', 10);
            closed_peaks(jj) = length(idx);
            
            
            this_z = zscore(out.wake{ii}(jj,:));
            [~, idx] = findpeaks(this_z, 'MinPeakHeight', 1.96, 'MinPeakDistance', 10);
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

exportgraphics(gcf, ['C:\Users\ecarm\Williams Lab Dropbox\Eric Carmichael\CIHR_2023_September\Assembly' filesep 'Assembly_Session_summary.pdf'], 'ContentType', 'vector')