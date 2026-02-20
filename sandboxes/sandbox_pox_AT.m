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


AT_mat = str2double(AT_tbl{:,5:17}); 
AT_sub = AT_tbl{:,1};  % subjects; 
AT_targets = AT_tbl.Properties.VariableNames(5:17); % target names
AT_depth = AT_tbl{:,4}; 
%% split the hippocampus into medial and lateral

midline = 2.5;


AT_mat_norm = AT_mat./AT_mat(:,8); 

s_list = unique(AT_sub); 

AT_mean = []; % make a new struct with the mean values 

for iS = 1:length(s_list)
    this_idx = contains(AT_tbl.Subject,s_list(iS));

    for ii = 1:length(AT_targets)

        if ismember(AT_targets{ii}, {'CA1', 'CA2', 'CA3', 'DG', 'DG'})
        
            this_idx = this_idx & (AT_depth < midline); 
            AT_mean.(AT_targets{ii})(iS) = mean(AT_mat_norm(this_idx, ii), 'omitmissing'); 

        else
            AT_mean.(AT_targets{ii})(iS) = mean(AT_mat_norm(this_idx, ii), 'omitmissing');

        end
    end
end