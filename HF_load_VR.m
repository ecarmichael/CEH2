function [vr] = HF_load_VR(fname)
%% HF_load_VR:  loads the .csv output from the Unity/arduino VR setup. 



%% newer version

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ";";

% Specify column names and types
opts.VariableNames = ["Time_s_", "ZPosition", "Event"];
opts.VariableTypes = ["double", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "Event", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "Event", "EmptyFieldRule", "auto");

vr_tbl = readtable(fname, opts); 
vr.info = []; 

vr.info.name = [fname(1:5) '_' fname(7:8)];
vr.info.date = strrep(fname(10:19), '-', '_');
vr.info.time = strrep(fname(21:28), '-', '_');


% position (continuous)
pos_idx = ismember(vr_tbl.Event, ""); 
vr.pos = tsd(vr_tbl.Time_s_(pos_idx), vr_tbl.ZPosition(pos_idx)'); 

% interpolate to get constant sampling rate
dt = 0.02; 
int_t = vr.pos.tvec(1):dt:vr.pos.tvec(end); 
vr.pos.data= interp1(vr.pos.tvec, vr.pos.data, int_t, 'linear', 'extrap');
vr.pos.tvec = int_t; 

% interval data
evts = unique(vr_tbl.Event(~ismember(vr_tbl.Event, ""))); 
vr.evt = ts; 

for iE  = length(evts):-1:1
    vr.evt.label{end+1} = evts{iE};
    vr.evt.t{end+1} = vr_tbl.Time_s_(ismember(vr_tbl.Event, evts{iE}));
end

% from the older version. 
% %% get the initialization information
% opts = delimitedTextImportOptions("NumVariables", 3, "Encoding", "UTF-8");
% % Specify range and delimiter
% opts.DataLines = [1, 22];
% opts.Delimiter = ",";
% opts.VariableNames = ["x0_05", "date", "x20250826104838"];
% opts.VariableTypes = ["double", "string", "string"];
% 
% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % Specify variable properties
% opts = setvaropts(opts, "date", "EmptyFieldRule", "auto");
% 
% vr_tbl = readtable(fname, opts); 
% vr.info = []; 
% % pull out each element
% d = vr_tbl{1,3};
% d = d{1}; 
% 
% vr.info.date = [d(1:4) '_' d(5:6) '_' d(7:8)];
% vr.info.time = [d(9:10) '_' d(11:12) '_' d(13:14)];
% 
% for ii = 3:length(vr_tbl.date)
%     vr.info.(vr_tbl.date{ii}) = vr_tbl{ismember(vr_tbl.date, vr_tbl.date{ii}),3};
%     vr.info.(vr_tbl.date{ii}) = vr.info.(vr_tbl.date{ii}){1};
%     if ~isnan(str2double(vr.info.(vr_tbl.date{ii})))
%         vr.info.(vr_tbl.date{ii}) = str2double(vr.info.(vr_tbl.date{ii}));
%     end
% end
% 
% 
% %% get the position informaiton
% 
% opts = delimitedTextImportOptions("NumVariables", 3, "Encoding", "UTF-8");
% opts.DataLines = [21 inf];
% opts.Delimiter = ",";
% opts.VariableNames = ["x0_05", "date","x20250826104838", "Var4", "Var5"];
% opts.VariableTypes = ["double", "string", "double", "double", "double"];
% 
% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % Specify variable properties
% opts = setvaropts(opts, "date", "EmptyFieldRule", "auto");
% 
% vr_tbl = readtable(fname, opts); 
% 
% 
% % position (continuous)
% pos_idx = ismember(vr_tbl.date, 'position'); 
% vr.pos = tsd(vr_tbl{pos_idx, {'x0_05'}}, vr_tbl{pos_idx, {'Var4'}}'); 
% 
% % interpolate to get constant sampling rate
% dt = 0.02; 
% int_t = vr.pos.tvec(1):dt:vr.pos.tvec(end); 
% vr.pos.data= interp1(vr.pos.tvec, vr.pos.data, int_t, 'linear', 'extrap');
% vr.pos.tvec = int_t; 
% 
% % interval data
% trig_idx = ismember(vr_tbl.date, 'trigger'); 
% vr.evt = ts({vr_tbl{trig_idx, {'x0_05'}}'}, {'trigger'}); 
