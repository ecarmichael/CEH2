function [vr] = HF_load_VR_TFC(fname)
%% HF_load_VR:  loads the .csv output from the Unity/arduino VR setup. 



%% newer version

opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["time", "encoderCount", "phase", "event", "context"];
opts.VariableTypes = ["double", "double", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["phase", "event", "context"], "EmptyFieldRule", "auto");

vr_tbl = readtable(fname, opts); 
vr.info = []; 

% vr.info.name = [fname(1:5) '_' fname(7:8)];
% vr.info.date = strrep(fname(10:19), '-', '_');
% vr.info.time = strrep(fname(21:28), '-', '_');

% set the time relative to the start of the protocol
vr_tbl.time = vr_tbl.time - vr_tbl.time(1);
vr_tbl.time = vr_tbl.time./1000; 

% position (continuous)
vr.pos = tsd(vr_tbl.time, abs(vr_tbl.encoderCount)'); 

% interpolate to get constant sampling rate
dt = 0.02; 
int_t = vr.pos.tvec(1):dt:vr.pos.tvec(end); 
[~, u_idx] = unique(vr.pos.tvec); 
vr.pos.data= interp1(vr.pos.tvec(u_idx), vr.pos.data(u_idx), int_t, 'linear', 'extrap');
vr.pos.tvec = int_t; 

%% interval data
vr.evt = ts; 

% phases
% phases = unique(vr_tbl.phase(ismember(vr_tbl.event, "-"))); 
% for iP = 1:length(phases)
%     vr.evt.label{end+1} = phases{iP};
%     p_idx = ismember(vr_tbl.phase, phases{iP}) & ismember(vr_tbl.event, "-"); 
%     vr.evt.t{end+1} = vr_tbl.time(p_idx);
% end

% other events
evts = unique(vr_tbl.event(~ismember(vr_tbl.event, ""))); 
for iE  = length(evts):-1:1
    vr.evt.label{end+1} = evts{iE};
    vr.evt.t{end+1} = vr_tbl.time(ismember(vr_tbl.event, evts{iE}));
end

% get the context states
Cons = unique(vr_tbl.context);
if isscalar(Cons)
    vr.evt.label{end+1} = Cons{iC};
    con_t = vr_tbl.time(ismember(vr_tbl.context, Cons(iC)));
    vr.evt.t{end+1} = [con_t(1) con_t(end)];
else
    for iC = 1:length(Cons)
        vr.evt.label{end+1} = [Cons{iC} '_on'];
        these_con = ismember(vr_tbl.context, Cons{iC});

        if iC == 1
            vr.evt.t{end+1} = [vr_tbl.time(1); vr_tbl.time(find(diff(these_con) ==-1)+1)];
        else
            vr.evt.t{end+1} = vr_tbl.time( find(diff(these_con) ==-1)+1);
        end


        vr.evt.label{end+1} = [Cons{iC} '_off'];
        if iC == length(Cons)
        vr.evt.t{end+1} = [vr_tbl.time(find(diff(these_con) ==1)+1);];

        else
        vr.evt.t{end+1} = vr_tbl.time(find(diff(these_con) ==1)+1);
        end

        if strcmp(vr_tbl.context(end), Cons{iC})
            vr.evt.t{end}(end+1) = vr_tbl.time(end); 
        end
    end
end


%% figure for varification
figure(1010)
clf

hold on

plot(vr.pos.tvec(1:end-1), MS_norm_range(diff(vr.pos.data(1,:)), -.5, .5))
c_ord = MS_linspecer(length(vr.evt.label));

for ii = 1:length(vr.evt.label)
    if strcmp(vr.evt.label{ii}, 'pos') || strcmp(vr.evt.label{ii}, 'Lick') 
        continue
    else
        % xline(vr.evt.t, 'color', c_ord(ii,:))
        % for jj = 1:length(vr.evt.t)
            % line(vr.evt.t{ii}, ones(size(vr.evt.t{ii}))*ii, 'color', c_ord(ii,:))
            plot([vr.evt.t{ii} vr.evt.t{ii}]', [(ones(size(vr.evt.t{ii}))*ii)-.5 (ones(size(vr.evt.t{ii}))*ii)+.5]', 'color', c_ord(ii,:))
            %  rectangle('position', [states.(s_list{ii}).tstart(jj), ii, (states.(s_list{ii}).tend(jj) - states.(s_list{ii}).tstart(jj)), 1], 'FaceColor', c_ord(ii,:));
        % end
    end
end

% s_list = replace(s_list, '_pos', ' +'); 
% s_list = replace(s_list, '_neg', ' -'); 

set(gca, 'YTick', 0:1:length(vr.evt.label), 'YTickLabel', ['mov' vr.evt.label])



%%  OLDER Set up the Import Options and import the data
% opts = delimitedTextImportOptions("NumVariables", 3);
% 
% % Specify range and delimiter
% opts.DataLines = [2, Inf];
% opts.Delimiter = ";";
% 
% % Specify column names and types
% opts.VariableNames = ["Time_s_", "ZPosition", "Event"];
% opts.VariableTypes = ["double", "double", "string"];
% 
% % Specify file level properties
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % Specify variable properties
% opts = setvaropts(opts, "Event", "WhitespaceRule", "preserve");
% opts = setvaropts(opts, "Event", "EmptyFieldRule", "auto");

% vr_tbl = readtable(fname, opts); 
% vr.info = []; 
% 
% vr.info.name = [fname(1:5) '_' fname(7:8)];
% vr.info.date = strrep(fname(10:19), '-', '_');
% vr.info.time = strrep(fname(21:28), '-', '_');
% 
% 
% % position (continuous)
% pos_idx = ismember(vr_tbl.Event, ""); 
% vr.pos = tsd(vr_tbl.Time_s_(pos_idx), vr_tbl.ZPosition(pos_idx)'); 
% 
% % interpolate to get constant sampling rate
% dt = 0.02; 
% int_t = vr.pos.tvec(1):dt:vr.pos.tvec(end); 
% vr.pos.data= interp1(vr.pos.tvec, vr.pos.data, int_t, 'linear', 'extrap');
% vr.pos.tvec = int_t; 
% 
% % interval data
% evts = unique(vr_tbl.Event(~ismember(vr_tbl.Event, ""))); 
% vr.evt = ts; 
% 
% for iE  = length(evts):-1:1
%     vr.evt.label{end+1} = evts{iE};
%     vr.evt.t{end+1} = vr_tbl.Time_s_(ismember(vr_tbl.Event, evts{iE}));
% end

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
