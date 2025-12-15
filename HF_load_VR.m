function [vr] = HF_load_VR(fname)
%% HF_load_VR:  loads the .csv output from the Unity/arduino VR setup. 




%% get the initialization information
opts = delimitedTextImportOptions("NumVariables", 3, "Encoding", "UTF-8");
% Specify range and delimiter
opts.DataLines = [1, 22];
opts.Delimiter = ",";
opts.VariableNames = ["x0_05", "date", "x20250826104838"];
opts.VariableTypes = ["double", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "date", "EmptyFieldRule", "auto");

vr_tbl = readtable(fname, opts); 
vr.info = []; 
% pull out each element
d = vr_tbl{1,3};
d = d{1}; 

vr.info.date = [d(1:4) '_' d(5:6) '_' d(7:8)];
vr.info.time = [d(9:10) '_' d(11:12) '_' d(13:14)];

for ii = 3:length(vr_tbl.date)
    vr.info.(vr_tbl.date{ii}) = vr_tbl{ismember(vr_tbl.date, vr_tbl.date{ii}),3};
    vr.info.(vr_tbl.date{ii}) = vr.info.(vr_tbl.date{ii}){1};
    if ~isnan(str2double(vr.info.(vr_tbl.date{ii})))
        vr.info.(vr_tbl.date{ii}) = str2double(vr.info.(vr_tbl.date{ii}));
    end
end


%% get the position informaiton

opts = delimitedTextImportOptions("NumVariables", 3, "Encoding", "UTF-8");
opts.DataLines = [21 inf];
opts.Delimiter = ",";
opts.VariableNames = ["x0_05", "date","x20250826104838", "Var4", "Var5"];
opts.VariableTypes = ["double", "string", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "date", "EmptyFieldRule", "auto");

vr_tbl = readtable(fname, opts); 


% position (continuous)
pos_idx = ismember(vr_tbl.date, 'position'); 
vr.pos = tsd(vr_tbl{pos_idx, {'x0_05'}}, vr_tbl{pos_idx, {'Var4'}}'); 

% interpolate to get constant sampling rate
dt = 0.02; 
int_t = vr.pos.tvec(1):dt:vr.pos.tvec(end); 
vr.pos.data= interp1(vr.pos.tvec, vr.pos.data, int_t, 'linear', 'extrap');
vr.pos.tvec = int_t; 

% interval data
trig_idx = ismember(vr_tbl.date, 'trigger'); 
vr.evt = ts({vr_tbl{trig_idx, {'x0_05'}}'}, {'trigger'}); 
