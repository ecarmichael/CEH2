%%sandbox POX AT level quantification


%% load the data

% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 21);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Subject", "Image_id", "Stain", "ML", "Vis", "CA2", "CA3", "DG", "CA1", "Rsc", "Sub", "STR", "PrSub", "PoSub", "MECSup", "MECDeep", "Var17", "Var18", "Var19", "ATMin", "ATMax"];
opts.VariableTypes = ["string", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var17", "Var18", "Var19"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Subject", "Image_id", "Var17", "Var18", "Var19"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "Stain", "TrimNonNumeric", true);
opts = setvaropts(opts, "Stain", "ThousandsSeparator", ",");

% Import the data
AT_tbl = readtable("/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/Histology/Pox_AT_levels - AT.csv", opts);


AT_mat = str2double(AT_tbl{1:19,5:17}); 
AT_sub = AT_tbl{1:19,1};  % subjects; 
AT_targets = AT_tbl.Properties.VariableNames(5:17); % target names
AT_depth = AT_tbl{1:19,4}; 
%% split the hippocampus into medial and lateral

midline = 2.5;


AT_mat_norm = AT_mat./AT_mat(:,8); 

s_list = unique(AT_sub); 

AT_mean = []; % make a new struct with the mean values 
AT_data_mat = []; 

for ii = 1:length(AT_targets)

    t_list{ii} = AT_targets{ii}; 

    for iS = 1:length(s_list)
        this_idx = contains(AT_sub,s_list(iS));

        if ismember(AT_targets{ii}, {'CA1', 'CA2', 'CA3', 'DG', 'DG'})

            this_idx = this_idx & (AT_depth < midline);
            AT_mean.(AT_targets{ii})(iS) = mean(AT_mat_norm(this_idx, ii), 'omitmissing');
            AT_data_mat(iS, ii) = mean(AT_mat_norm(this_idx, ii), 'omitmissing');

        else
            AT_mean.(AT_targets{ii})(iS) = mean(AT_mat_norm(this_idx, ii), 'omitmissing');
            AT_data_mat(iS, ii) = mean(AT_mat_norm(this_idx, ii), 'omitmissing');
        end
    end
end

% reorder
AT_data = []; 

AT_data(1,:) = AT_mean.STR; 
AT_data(2,:) = AT_mean.Sub; 
AT_data(3,:) = AT_mean.PrSub; 
AT_data(4,:) = AT_mean.PoSub; 
AT_data(5,:) = AT_mean.MECSup; 
AT_data(6,:) = AT_mean.MECDeep;
AT_data(7,:) = AT_mean.CA1; 
AT_data(8,:) = AT_mean.CA2; 
AT_data(9,:) = AT_mean.CA3; 
AT_data(10,:) = AT_mean.DG;
AT_data(11,:) = AT_mean.Rsc; 
AT_data(12,:) = AT_mean.Vis;

AT_data  = AT_data'; 

AT_targets_sort  = {'STR', 'Sub', 'PreSub', 'PostSub', 'MEC_{Sup}', 'MEC_{Deep}', 'CA1', 'CA2', 'CA3', 'DG', 'RSC', 'Vis'}; 

AT_tbl_out  = cell2table(num2cell(AT_data), 'VariableNames',AT_targets_sort, 'RowNames',s_list); 

%% plot points for each subject and region. 

% t_list = fieldnames(AT_mean); 

% AT_targets_sort  = {'STR', 'MECSup', 'MECDeep', 'Sub', 'PrSub', 'PostSub', 'RSC', 'Vis', 'CA1', 'CA2', 'CA3', 'DG'}; 

c_ord = [1 1 1; spring(5); winter(4); prism(2)]; 

figure(10)
clf

hold on
for ii = 1:length(AT_targets_sort)

    this_data = AT_data(:,ii); 

    x_vals = repmat(ii, length(this_data),1); 
    % x_vals = x_vals+ sort(MS_randn_range(length(x_vals), 1, -.1, .1));

    scatter(x_vals, this_data*100,100,  'filled', 'MarkerFaceColor',c_ord(ii,:))

    bar(ii, mean(this_data*100, 'omitmissing'), 'facecolor', 'none')

end

b_ord = nebula(size(AT_data,1)); 

x_val = 1:length(AT_data); 
for ii = 1:size(AT_data,1)

    nan_idx = isnan(AT_data(ii,:)); 

    plot(x_val(~nan_idx), AT_data(ii,~nan_idx)*100, '.--', 'color', b_ord(ii,:))

end
ylabel({'Mean AT fluoresence', 'realtive to STR'})
set(gca, 'XTick',1:length(AT_targets_sort), 'XTickLabel', AT_targets_sort)

yline(100)

writetable(AT_tbl_out, ["/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/Histology/AT_level_tbl.xls"])


%% same thing but relative to Sub

AT_data_r_sub = AT_data./AT_data(:,2); 

figure(20)
clf

hold on
for ii = 1:length(AT_targets_sort)

    this_data = AT_data_r_sub(:,ii); 

    x_vals = repmat(ii, length(this_data),1); 
    % x_vals = x_vals+ sort(MS_randn_range(length(x_vals), 1, -.1, .1));

    scatter(x_vals, this_data*100,100,  'filled', 'MarkerFaceColor',c_ord(ii,:))

    bar(ii, mean(this_data*100, 'omitmissing'), 'facecolor', 'none')

end

b_ord = nebula(size(AT_data_r_sub,1)); 

x_val = 1:length(AT_data_r_sub); 
for ii = 1:size(AT_data_r_sub,1)

    nan_idx = isnan(AT_data_r_sub(ii,:)); 

    plot(x_val(~nan_idx), AT_data_r_sub(ii,~nan_idx)*100, '.--', 'color', b_ord(ii,:))

end
ylabel({'Mean AT fluoresence', 'realtive to STR'})
set(gca, 'XTick',1:length(AT_targets_sort), 'XTickLabel', AT_targets_sort)

yline(100)
