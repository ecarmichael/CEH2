function [out] = MS_NOL_quick_score(fname, scale_f)

if nargin < 1
    scale_f = 1; 
    d = dir('*markers.csv');
    fname = d.name; 
end

if isempty(fname)
        d = dir('*markers.csv');
    fname = d.name; 
end
%% quick Chronotate score

opts = delimitedTextImportOptions("NumVariables", 4);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["marker_time", "marker_type", "object_id", "marker_frame"];
opts.VariableTypes = ["double", "string", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "marker_type", "EmptyFieldRule", "auto");

% Import the data
tbl = readtable(fname, opts);

%% cound the start and tops. 

tbl.marker_time = tbl.marker_time *scale_f; 

    % get all the object epochs
    obj_1_idx = find((tbl.object_id == 1 ) & contains(tbl.marker_type, 'start')); 
    obj_2_idx = find((tbl.object_id == 2 ) & contains(tbl.marker_type, 'start')); 

    % remove any blocks starting after 600seconds
    rm_idx = tbl.marker_time(obj_1_idx) > 300; 
    obj_1_idx(rm_idx) = []; 
    
    rm_idx = tbl.marker_time(obj_2_idx) > 300; 
    obj_2_idx(rm_idx) = []; 

    % loop over events and get the duration
    obj_1 = []; 
    for jj = length(obj_1_idx):-1:1

        if  tbl.marker_time(obj_1_idx(jj)+1) > 300
            obj_1(jj) = 300 - tbl.marker_time(obj_1_idx(jj));
        else
            obj_1(jj) = tbl.marker_time(obj_1_idx(jj)+1) - tbl.marker_time(obj_1_idx(jj)); 
        end
    end

    obj_2 = []; 
    for jj = length(obj_2_idx):-1:1
        obj_2(jj) = tbl.marker_time(obj_2_idx(jj)+1) - tbl.marker_time(obj_2_idx(jj)); 
    end


    data_out.Obj_1_n = length(obj_1); 
    data_out.Obj_1_t = sum(obj_1); 
    data_out.Obj_2_n = length(obj_2); 
    data_out.Obj_2_t = sum(obj_2); 
    data_out.ratio_t = (sum(obj_2) - sum(obj_1))/(sum(obj_2)+sum(obj_1)); 
    data_out.ratio_n = (length(obj_2) - length(obj_1))/(length(obj_2)+length(obj_1)); 

    fprintf("Obj 1 n: %.0f  t:%.2f    |  Obj 2 n: %.0f  t:%.2f   | <strong>%.1f</strong> <strong>(%.1f)</strong>\n", length(obj_1), sum(obj_1), length(obj_2), sum(obj_2), data_out.ratio_t, data_out.ratio_n)
