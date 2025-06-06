function [csc_out, ms_out]  = CE_sandbox_load_NLX(nlx_dir, CA,chan_id evt_idx)


%% init
if nargin < 4

    evt_idx = []; 
    fprintf('<strong>%s</strong>: no ''evt_idx'' specified using closest match', mfilename)
end


if isstring(CA)
    CA_dir = CA; 
        fprintf('<strong>%s</strong>: ''CA'' is a string, using it as a directory', mfilename)

elseif isstruct(CA)
    CA_type = 'ms';
        fprintf('<strong>%s</strong>: ''CA'' input is a struct aligning to pre and post', mfilename)
    ms_pre = CA.pre.ms; 
    ms_post = CA.post.ms; 

end


% nlx_dir = 'K:\NLX_CA test\2023-09-29_TTL_Test'; 
% 
% CA_dir = {'K:\NLX_CA test\12_35_04_1827\My_V4_Miniscope',...
%     'K:\NLX_CA test\12_35_36_1828\My_V4_Miniscope'}; 
%% grab some meta_data
cd(nlx_dir); 

cfg_csc = [];
cfg_csc.fc = {chan_id}; 
csc = MS_LoadCSC(cfg_csc);


% parts = strsplit(cd, filesep);
% parts = strsplit(parts{end}, '_'); 
% meta = [];
% meta.proj = parts{1};
% meta.sub1 = parts{2};
% meta.sub2 = parts{3};
% meta.sub3 = parts{4};
% meta.date = parts{5};
% meta.rec = parts{6};


evts = LoadEvents([]);

% split into the two longest recording blocks

rec_len = []; 
start_idx = find(contains(evts.label, 'Starting Recording'));

for ii = length(evts.t{start_idx}):-1:1

    rec_len(ii) = evts.t{start_idx+1}(ii) - evts.t{start_idx}(ii); 

end

if length(rec_len) > 2
    error('More than 2 recording blocks (%f), fix this code if true', length(rec_len)); 
    % [len, idx] = sort(rec_len); 
end

pre_evts = restrict(evts,evts.t{start_idx}(1), evts.t{start_idx+1}(1)); 
post_evts = restrict(evts,evts.t{start_idx}(2), evts.t{start_idx+1}(2)); 

rm_idx = contains(pre_evts.label, 'St');
pre_evts.t(rm_idx) = [];
pre_evts.label(rm_idx) = [];

rm_idx = contains(post_evts.label, 'St');
post_evts.t(rm_idx) = [];
post_evts.label(rm_idx) = [];
%% compare the number of timestamps between NLX and CA

if ~isempty(evt_idx)
    % ms_pre
    csc_out.pre = []; 
    csc_out.pre.evts = pre_evts;
    csc_out.pre.evts.t =[];
    csc_out.pre.evts.t{1} =sort(unique([pre_evts.t{1}, pre_evts.t{2}]));
    csc_out.pre.evts.label = [];
    csc_out.pre.evts.label = 'ms_ttl';
    % restrict the csc to the first and last TTl
    csc_out.pre.csc = restrict(csc,csc_out.pre.evts.t{1}(1), csc_out.pre.evts.t{1}(end));   
    %csc_out.pre.csc.tvec = csc_out.pre.csc.tvec - csc_out.pre.csc.tvec(1); 
    %csc_out.pre.evts.t{1} = csc_out.pre.evts.t{1} - csc_out.pre.evts.t{1}(1)-evts.t{start_idx}(1); 



    csc_out.post = []; 
    csc_out.post.evts = post_evts;
    csc_out.post.evts.t =[];
    csc_out.post.evts.t{1} =sort(unique([post_evts.t{1}, post_evts.t{2}]));
    csc_out.post.evts.label = [];
    csc_out.post.evts.label = 'ms_ttl';
    % restrict the csc to the first and last TTl
    csc_out.post.csc = restrict(csc,csc_out.post.evts.t{1}(1), csc_out.post.evts.t{1}(end));  

  % fprintf('<strong>%s</strong>: NLX = %.0f  | MS = %.0f (<strong>%.0f</strong>)  time diff: %.2f\n',...
        % mfilename, length(csc_out.pre.evts.t{1}), length(ms_pre.time),length(csc_out.pre.evts.t{1})-length(ms_pre.time),...
        % (csc_out.pre.evts.t{1}(end) - csc_out.pre.evts.t{1}(1)) - 




else
    
    fprintf('PRE\n')
 for ii = 1:size(pre_evts.t,2)
    fprintf('<strong>NLX evt %d</strong>: NLX = %.0f  | MS = %.0f (<strong>%.0f</strong>) \n',ii, length(pre_evts.t{ii}), length(ms_pre.time),length(pre_evts.t{ii})-length(ms_pre.time))
 end

    % make a matrix of all of the on and off times
    n_ttl = []; 
    for ii = size(pre_evts.t,2):-1:1
            for jj = size(pre_evts.t,2):-1:1
                n_ttl(ii, jj) = length(sort(unique([pre_evts.t{ii};pre_evts.t{jj}])) );

            end
    end
    n_ttl = n_ttl - length(ms_pre.time)





end





%%

rm_idx = contains(evts.label, 'St');
evts.t(rm_idx) = [];
evts.label(rm_idx) = [];

count = 1; TTLs = [];
% for ii = 1:2:length(evts.t)
    TTLs{1} = sort(unique([evts.t{1} evts.t{2}]));     
    TTLs{2} = sort(unique([evts.t{3} evts.t{4}])); 
    TTLs{3} = sort(unique([evts.t{5} evts.t{6}])); 

%     count = count+1;
% end




%% load the timestamps from the miniscope. 
% json_f = dir('*.json');
% 
% fid = fopen(json_f.name);
%     raw = fread(fid,inf);
%     str = char(raw');
%     Exp_json = jsondecode(str); 
%     fclose(fid);
fnum = []; tvec = [];
for ii = 1:length(CA_dir)
    cd(CA_dir{ii})
    TS = readtable('timeStamps.csv');
    fnum{ii} = TS.FrameNumber;
    tvec{ii} = TS.TimeStamp_ms_ ./1000; % convert to seconds
    
end
    
%% compare the evts and the TS

for ii = 1:length(tvec)
        fprintf('TS %0.0f = %0.0f fs = %0.2f\n',ii, length(tvec{ii}), 1/mode(diff(tvec{ii})))

end

for ii = 1:length(TTLs)
            fprintf('TTLs %0.0f = %0.0f fs = %0.2f\n',ii, length(TTLs{ii}), 1/mode(diff(TTLs{ii})))
end

% for ii = 1:length(evts.t)
%             fprintf('evts %0.0f = %0.0f fs = %0.2f\n',ii, length(evts.t{ii}), 1/mode(diff(evts.t{ii})))
% end

%% test plot

figure(101)
clf

hold on
plot(csc.tvec, csc.data)

y_lim = ylim; 

c_ord = linspecer(length(evts.t));
    offset = y_lim(1) + min(csc.data); 

for ii = 1:length(TTLs)
   plot([TTLs{ii}; TTLs{ii}], [ones(1,length(TTLs{ii}))*(offset*ii); ones(1,length(TTLs{ii}))*(offset*ii)+offset],'-', 'color', c_ord(ii,:))
%        plot([evts.t{ii}; evts.t{ii}], [ones(1,length(evts.t{ii}))*(offset*ii); ones(1,length(evts.t{ii}))*(offset*ii)+offset],'-', 'color', c_ord(ii,:))

end