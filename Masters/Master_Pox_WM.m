%% Master collector script for Pox WM data




%% load the learning rates for each cohort. 



all.Jan25 = readtable('POX_WM - Jan25_Cohort.csv'); 





%% loop over cohorts and grab all the subjects
c_ord = MS_linspecer(4); 

coh = fieldnames(all); 

for iC = 1:length(coh)

    subs = all.(coh{iC}).ID; 
    keep_idx = contains(subs, 'Pox'); 

    this_tab = all.(coh{iC}); 
    this_tab(~keep_idx,:) = []; 

    this_data = table2array(this_tab(:,4:23));
    this_day_mean = []; 
    for ii = size(this_data,1):-1:1
            this_day_mean(ii,1:5) = [mean(this_data(ii,1:4)) mean(this_data(ii,5:8)), mean(this_data(ii,9:12)), mean(this_data(ii,13:16)), mean(this_data(ii,17:20))];
    end



    data.(coh{iC}).sub = table2cell(this_tab(:,1)); 

    data.(coh{iC}).geno= logical(table2array(this_tab(:,2))); 

    data.(coh{iC}).data =  this_data; 

    data.(coh{iC}).day_mean =  this_day_mean; 

    data.(coh{iC}).day_mean_ctrl =  this_day_mean(data.(coh{iC}).geno == 0,:); 

    data.(coh{iC}).day_mean_tau =  this_day_mean(data.(coh{iC}).geno == 1,:);

    data.(coh{iC}).cross = table2array(this_tab(:,24)); 




    figure(iC)
    clf
    title(coh{iC})
    subplot(1, 4, 1:3)
    cla
    hold on

    % for ii = 1:length(data.(coh{iC}).sub)
    %     plot(1:5, this_day_mean(ii,:))
    % 
    % end
    plot(1:5, mean(data.(coh{iC}).day_mean_ctrl), 'color', c_ord(1,:), 'LineWidth',2)
    plot(1:5, mean(data.(coh{iC}).day_mean_tau), 'color', c_ord(2,:), 'LineWidth',2)

    plot(1:5, mean(this_day_mean), 'k', 'LineWidth',2)
    % leg = data.(coh{iC}).sub; 
    % leg{end+1} = 'mean';
    legend({'Ctrl', 'Tau', 'All'}, 'Box','off')
    ylabel('Mean escape latency')

    set(gca, 'XTickLabel', {'Day 1' 'Day2 ' 'Day 3' 'Day 4' 'Day 5'}, 'XTickLabelRotation', 45)

    subplot(1,4,4)
        
    b = MS_bar_w_err(data.(coh{iC}).cross(data.(coh{iC}).geno == 0),data.(coh{iC}).cross(data.(coh{iC}).geno == 1),c_ord(1:2,:), 1, 'ttest2' );
    b(1).FaceColor = 'none';
    b(1).EdgeColor = 'k';

    scatter(1+sort(MS_randn_range(length(data.(coh{iC}).cross(data.(coh{iC}).geno == 0)), 1, -.1, .1)), data.(coh{iC}).cross(data.(coh{iC}).geno == 0),25,  c_ord(1,:), 'filled')
    scatter(2+sort(MS_randn_range(length(data.(coh{iC}).cross(data.(coh{iC}).geno == 1)), 1, -.1, .1)), data.(coh{iC}).cross(data.(coh{iC}).geno == 1),25,  c_ord(2,:), 'filled')


end

%% simple plots

plot(this_tab(2:end,:))

