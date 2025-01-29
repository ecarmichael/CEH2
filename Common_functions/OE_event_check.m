%% OE_event check


%% load some data


evts = OE_LoadEvents();

cfg_csc = [];
cfg_csc.fc =  {'CH1'};%cfg.csc_chan;
% cfg_csc.desired_sampling_frequency = 2000;
csc = OE_old_csc2TSD(cfg_csc);

%% remove jitter
for ii = 1:length(evts.label)
    
%    temp_i = evts.t{ii};
%    
%    f_p = find(diff(temp_i)>.5)+1; % only get events with more than one 10ms between them 
%     
%    evts.t{ii} = evts.t{ii}([1 f_p']);


this_iv = iv;
this_iv.tstart = evts.t{ii}(1:2:end); 
this_iv.tend = evts.t{ii}(2:2:end); 

this_iv.usr.len = this_iv.tend - this_iv.tstart; 



this_iv = SelectIV([], this_iv, this_iv.usr.len > .1); 

out = []; 
for jj = 1:length(this_iv.tstart)
    out = [out this_iv.tstart(jj) this_iv.tend(jj)]; 
end


evts.t{ii} = out; 

end

%% plot csc and events

c_ord =  MS_linspecer(length(evts.label)+1);


figure(1914)
clf

hold on
plot(csc.tvec, csc.data)
% 
y_l = ylim; 
h = y_l(2) - y_l(1);  

for ii = [1 3]
%     
%        
       for jj = 1:2:length(evts.t{ii})
           rectangle('position', [evts.t{ii}(jj),y_l(1), evts.t{ii}(jj+1) - evts.t{ii}(jj), h], 'FaceColor', [c_ord(ii,:) .4])
       end
end

%% get the start times realtive to the start of the sleep phase

% convert to seconds

t_l = evts.t{1}; 
for ii = 1:length(evts.t{1})

    t = evts.t{1}(ii);
    
    t_min = fix(t/ 60); 
    t_s = mod(abs(t/60),1)*60; 

    t_m(ii) = t_min + t_s/100;
end

% write to a file
parts = strsplit(cd, filesep);
if ~isunix
parent_path = strjoin(parts(1:end-1), [filesep filesep]);
else
    parent_path = strjoin(parts(1:end-1), [filesep ]);
end

fid = fopen([parent_path filesep 'light_times.txt'], 'w');



c = 0; 
for ii  = 1:2:length(t_l)
    c=  c+1;
    fprintf(fid,'%.0f:  %.2fmin - %.1fs duration\n', c, t_m(ii), (t_l(ii+1) - t_l(ii))*60);
        fprintf('%.0f:  %.2fmin - %.1fs duration\n', c, t_m(ii), (t_l(ii+1) - t_l(ii))*60);

end


fclose(fid);







