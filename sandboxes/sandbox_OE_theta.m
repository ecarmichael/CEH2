%% sandbox_PV_opto check



%% load the OE csc for specific animals

parts = strsplit(pwd, filesep);

fname = parts{end};
s_idx = strfind(fname, '2023');

sub = fname(1:3); 
sess = strrep(fname(s_idx:s_idx+9), '-', '_'); 

if strcmpi(sub, 'PV1')
    LFP_chan = 'CH63';
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

keep_t =evts.t{1}(keep_idx); 
LED_IV = iv(keep_t(1:2:end), keep_t(2:2:end));

% LED_pace = iv(evts.t{1}(1)*0:20:200-10, evts.t{1}(1):10:200); 

%%

figure(101)
clf
hold on
for ii = 1:length(LED_IV.tstart)
    rectangle('Position',[LED_IV.tstart(ii), min(csc.data), 10, max(csc.data)-min(csc.data)], 'FaceColor', 'r')
end
plot(csc)


%%
cfg_psd = []; 
cfg_psd.hann_win = 2^11; 

all_ppx_LED = [];
for ii =length(LED_IV.tstart):-1:1
   this_data = restrict(csc, LED_IV.tstart(ii), LED_IV.tend(ii)); 
   
   [ppx, f] = MS_get_psd(cfg_psd, this_data);
    
   all_ppx_LED(:,ii) = ppx; 
       all_BP_LED(ii) = 10*log10(bandpower(ppx, f, [6 10], 'psd')); 
     all_BP_norm_LED(ii) = 10*log10(bandpower(ppx, f, [1 55], 'psd')); 

end

all_LED_m = mean(all_ppx_LED,2); 

%% same for the laser off
all_ppx_off = [];
for ii = length(LED_IV.tstart):-1:1
   this_data = restrict(csc, LED_IV.tstart(ii)+10, LED_IV.tend(ii)+10); 
   
   [ppx, f] = MS_get_psd(cfg_psd, this_data);
    
   all_ppx_off(:,ii) = ppx; 
   all_BP_off(ii) = 10*log10(bandpower(ppx, f, [6 10], 'psd')); 

     all_BP_norm_off(ii) = 10*log10(bandpower(ppx, f, [1 55], 'psd')); 
    
end
all_off_m = mean(all_ppx_off,2); 



%% plot

figure(101)
clf
subplot(2,2,1)
cla
hold on
plot(f, 10*log10(all_LED_m), 'color', [0 .8 0.2])
plot(f, 10*log10(all_off_m), 'color', [.8 .8 .8])

xlim([0 100])



subplot(2,2,3)
cla
b = bar(1:2,[mean(all_BP_LED)/mean(all_BP_norm_LED); mean(all_BP_off)/mean(all_BP_norm_off)]); 
b.FaceColor = 'flat';
b.CData(1,:) = [0 .8 0.2];
b.CData(2,:) = [0.8 .8 0.8];


%% check the trajectory

up_dir = fullfile(cd, '..')
cd(up_dir)

%% move up to the miniscope dir. 
mini_dir = dir(fullfile(cd, '*minicam*'))
cd([mini_dir.folder filesep mini_dir.name])


%%
pos = MS_DLC2TSD(cd);


figure(101)
subplot(2,2,[2 4])
cla
plot(pos.data(1,:), pos.data(2,:), '.')