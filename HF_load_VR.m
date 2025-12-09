function [VR_out] = HF_load_VR(fname)
%% HF_load_VR:  loads the .csv output from the Unity/arduino VR setup. 




%% get the initialization information
opts = delimitedTextImportOptions("NumVariables", 3, "Encoding", "UTF-8");
% Specify range and delimiter
opts.DataLines = [1, 22];
opts.Delimiter = ",";
opts.VariableNames = ["x0_05", "date", "x20250826104838"];
opts.VariableTypes = ["double", "categorical", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "date", "EmptyFieldRule", "auto");

vr_tbl = readtable(fname, opts); 

% pull out each element
d = vr_tbl{1,3};
d = d{1}; 

vr.info.date = [d(1:4) '_' d(5:6) '_' d(7:8)];
vr.info.time = [d(1:4) '_' d(5:6) '_' d(7:8)];