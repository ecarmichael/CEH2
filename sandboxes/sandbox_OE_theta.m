%% sandbox_PV_opto check



%% load the OE csc for specific animals

parts = strsplit(pwd, filesep);

fname = parts{end};
s_idx = strfind(fname, '2023');

sub = fname(1:3); 
sess = strrep(fname(s_idx:s_idx+9), '-', '_'); 

if strcmpi(sub, 'PV1')
    LFP_chan = 'CH61';
elseif strcmpi(sub, 'PV2')
    LFP_chan = 'CH61';



end

cfg.TTL_frame = '6'; 
cfg.TTL_trig = '1';

cd('Record Node 106/')
%% load the events

evts = OE_LoadEvents; 

if sum(contains(evts.label, cfg.TTL_trig)) == 1
    evts.label{(contains(evts.label, cfg.TTL_frame))} = 'Rec_TTL'; 
    evts.label{(contains(evts.label, cfg.TTL_trig))} = 'LED_TTL'; 
    
    for ii = 1:length(evts.label)
        dur(ii) = evts.t{ii}(end) - evts.t{ii}(1); 
    end
    
    [~, max_TTL] = max(dur); 
    start_t = evts.t{max_TTL}(1);
    end_t = evts.t{max_TTL}(end); 
    
else
    error('No trigger TTL found. Required for synchronizing data')
end

fprintf('<strong>%s</strong>: OE frame TTL (<strong>''%s''</strong>), duration = %.0fsec (%2.1fhr)\n', mfilename,evts.label{max_TTL}, end_t - start_t, (end_t - start_t)/60/60)



%% cload the CSC


cfg_csc = [];
cfg_csc.fc =  {LFP_chan};%cfg.csc_chan;
cfg_csc.desired_sampling_frequency = 2000;
csc = OE_old_csc2TSD(cfg_csc);

%% extract the LED on IVs
keep_idx = diff(evts.t{1}) > 8;
keep_idx = [1 ; find(keep_idx(1:end))+1]';

LED_IV = iv(evts.t{1}(keep_idx), evts.t{1}(keep_idx)+10);

LED_pace = evts.t{1}(1)*0:10:200; 

%%

figure(101)
clf
hold on
for ii = 1:length(LED_IV.tstart)
    rectangle('Position',[LED_IV.tstart(ii), min(csc.data), 10, max(csc.data)-min(csc.data)], 'FaceColor', 'r')
end
plot(csc)


%%


MS_detect_SWR_JC()