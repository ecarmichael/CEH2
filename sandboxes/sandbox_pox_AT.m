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



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%    NEW BATCH   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

opts = delimitedTextImportOptions("NumVariables", 7);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Image", "Annotation", "Mean", "Median", "Min", "Max", "PixelCount"];
opts.VariableTypes = ["categorical", "string", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Annotation", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Image", "Annotation"], "EmptyFieldRule", "auto");



data_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/Histology/Pox_batch/project/results'; 

f_names = dir([data_dir filesep '*.csv']); 

all_tbl = []; 

for ii = 1:length(f_names)

    % load the data 
    this_tbl = readtable(f_names(ii).name); 

    % restrict the name to the mouse ID. 
    for jj = length(this_tbl.Image):-1:1

        % rename
        d_idx = strfind(this_tbl.Image{jj}, '2026')-1;

        this_tbl.slide{jj} = this_tbl.Image{jj}(end-6:end);

        this_tbl.Image{jj} = this_tbl.Image{jj}(1:d_idx); 
    end

% normalize to the str
    str_idx = find(contains(upper(this_tbl.Annotation), 'STR')); 

    this_norm_tbl = this_tbl; 

    % for jj = length(this_tbl.Image):-1:1
        this_tbl.Mean_norm = this_tbl.Mean./this_norm_tbl.Mean(str_idx);

    % end

    all_tbl = [all_tbl;  this_tbl]; 
end

all_tbl = movevars(all_tbl, 'slide', 'After', 'Image'); 
all_tbl = movevars(all_tbl, 'Mean_norm', 'After', 'Mean'); 
all_tbl.Annotation = lower(all_tbl.Annotation); 

%% get the means for each animals

sub_list = unique(all_tbl.Image);
reg_list = unique(all_tbl.Annotation);

data_out = NaN(length(sub_list), length(reg_list)); 

for ii = 1:length(sub_list)

s_idx = contains(all_tbl.Image, sub_list{ii}); 

    for rr = 1:length(reg_list)

        r_idx = contains(all_tbl.Annotation, reg_list{rr}); 

        this_idx =s_idx & r_idx; 

        data_out(ii, rr) = mean(all_tbl.Mean_norm(this_idx), 'omitmissing'); 

    end
end

keep_idx = ismember(lower(reg_list), {'str', 'ca1', 'sub'}); 

data_out(:, ~keep_idx) = []; 

ctrl_idx = logical([1 1 0 0 0 0 1 0 1]); 

trace = [.62, 1, .10, .4, .88, .95, .53, .67, .52];


c_ord = MS_linspecer(length(ctrl_idx)); 

figure(201)
clf

subplot(2,4,[1:2 5 6])
hold on
for ii = 1:size(data_out, 1)

    if ctrl_idx(ii)
        scatter(1:size(data_out,2), data_out(ii,:), 100,[.7 .7 .7],'s',  "filled")
        plot(1:size(data_out,2), data_out(ii,:), 'color', [.7 .7 .7])
    else
        scatter(1:size(data_out,2), data_out(ii,:), 100,c_ord(ii,:),'s',  "filled")
        plot(1:size(data_out,2), data_out(ii,:), 'color', c_ord(ii,:))
    end
end

set(gca, 'xtick', 1:size(data_out,2), 'XTickLabel', reg_list(keep_idx))
xlim([0 size(data_out,2)+1])


subplot(2, 4, 3)
[hb, eb, sc, p, stats] = MS_bar_w_err(trace(ctrl_idx), trace(~ctrl_idx),[.7 .7 .7; c_ord(1,:)], 1, 'ttest2'); 
set(gca, 'color', 'w')


subplot(2, 4, 4)
sub_idx = find(ismember(reg_list(keep_idx), 'sub')); 
[hb, eb, sc, p, stats] = MS_bar_w_err(data_out(ctrl_idx,sub_idx), data_out(~ctrl_idx,sub_idx),[.7 .7 .7; c_ord(1,:)], 1, 'ttest2'); 
set(gca, 'color', 'w')

subplot(2, 3, 6)
hold on

scatter(data_out(~ctrl_idx,sub_idx), trace(~ctrl_idx), 50, c_ord(1,:), 'filled')
scatter( data_out(ctrl_idx,sub_idx), trace(ctrl_idx), 50, [.7 .7 .7], 'filled')

ylabel('Trace Freezing %')
xlabel('Sub pTau')

set(gca, 'XScale', 'log')


