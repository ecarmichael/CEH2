function evts = MS_LoadEvents
%% MS_LoadEvents: simplified version of LoadEvents from the vandermeerlab codebase. 
%
%
%
%    Inputs: 
%    - 
%
%
%
%    Outputs: 
%    - evts: [struct] contains all of the events and time stamps in the ts
%    format. 
%
%
%
%
% EC 2023-10-30   this simple version 
%  MvdM 2014-06-18, original 
%
%
%% initialize

e_dir = dir('Events.nev');

if isempty(e_dir)
    error('''no Events.nev'' file found')
end



%% load using NLX mex file. 



[EVTimeStamps, EventIDs, ~, EVExtras, EventStrings, EVHeader] = Nlx2MatEV(e_dir.name,[1 1 1 1 1],1,1,[]);

%% pull out the unique events

eventList = unique(EventStrings);


eventList(contains(eventList, 'CSC')) = []; % remove extra CSC Reference change events.  Odd case. 
%% for each event type get the times

evts= ts();
% evts.type = 'ts';
% evts.t = [];
% evts.label = [];
% evts.cfg.hdr = [];

for iE = 1:length(eventList)
   
    ev_string = eventList{iE};

    ev_id = strncmp(eventList{iE},EventStrings,length(eventList{iE}));
    ev_t = EVTimeStamps(ev_id)*10^-6;
    
    % check if this eventLabel already exists, if so append (not create new)
    label_idx = strmatch(EVHeader{iE},evts.label,'exact');
    
    if isempty(label_idx) % label doesn't exist, create new
        
        evts.t{iE} = ev_t;
        evts.label{iE} = eventList{iE};
        
    else % label exists, append and sort (sort is helpful for later use with nearest_idx3)
        
        evts.t{label_idx} = sort(cat(2,evts.t{label_idx},ev_t));
                
    end
    
end

%% conver the header
evts.cfg.hdr = [];

for hline = 1:length(EVHeader)
    
    line = strtrim(EVHeader{hline});
    
    if isempty(line) || ~strcmp(line(1),'-') % not an informative line, skip
        continue;
    end
    
    % if we are reading an older version header with colons, remove them
    colon_idx = strfind(line,':');
    if ~isempty(colon_idx)
       line(colon_idx) = []; 
    end
    
    % This expression captures the first chunk of alphanumeric characters and puts it into
    % <key> then puts whatever is to the right of it into <val>. If there is only one
    % character chunk (e.g., missing value) then it returns <val> as empty.
    a = regexp(line(2:end),'(?<key>^\S+)\s+(?<val>.*)|(?<key>\S+)','names');

    % deal with characters not allowed by MATLAB struct
    if strcmp(a.key,'DspFilterDelay_ï¿½s') | strcmp(a.key,'DspFilterDelay_µs')
        a.key = 'DspFilterDelay_us';
    end
    
    evts.cfg.hdr = setfield(evts.cfg.hdr,a.key,a.val);
    
    % convert to double if possible
    if ~isnan(str2double(a.val))
        evts.cfg.hdr = setfield(evts.cfg.hdr,a.key,str2double(a.val));
    end
    
end


