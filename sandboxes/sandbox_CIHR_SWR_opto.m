%% sandbox_CIHR_SWR_opto


evt_bin_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/H2b32026-03-15_00-05-23_dSub_SWR/Record Node 112/experiment1/recording1/events/Intan_RHD_USB-158.Rhythm Data/TTL';

csc_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/H2b32026-03-15_00-05-23_dSub_SWR/Record Node 117'; 
swr_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/H2b32026-03-15_00-05-23_dSub_SWR/Record Node 143';
bin_dir =     '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/H2b32026-03-15_00-05-23_dSub_SWR/Record Node 112/experiment1/recording1/continuous/Intan_RHD_USB-158.Rhythm Data';

% load the binary files


% evts = OE_load_binary_evts(evts_dir, );
ts_prime = readNPY([bin_dir filesep 'timestamps.npy']); 

% cd(evt_bin_dir)
evts = OE_load_binary_evts(evt_bin_dir, ts_prime(1));

cd(swr_dir)
evts_csc = OE_LoadEvents(); 


cd(csc_dir)
cfg.fc = {'CH11', 'CH12', 'CH14'};
csc = MS_LoadCSC_OE(cfg);

csc.tvec = csc.tvec - ts_prime(1); 


for ii = 1:length(evts_csc)

evts_csc.t{ii} = evts_csc.t{ii} - ts_prime(1); 
end

%% swr detection

csc.data(1,:)= mean(csc.data,1); 
csc.label{1} = 'swr_mean'; 

swr = MS_SWR_detector(csc, csc.label{1}); 

%% plot to see if everything is aligned
c_ord = MS_linspecer(3); 

figure(101)
clf
hold on

% plot(csc.tvec, csc.data(1,:),'k')
cfg_plot.target = csc.label{1};   
PlotTSDfromIV(cfg_plot, swr, csc)

vline(evts.t{2}(1,:), 'r');
vline(evts.t{2}(2,:), 'm');
vline(evts.t{3}(1,:), 'b');
vline(evts.t{3}(2,:), 'c');

vline(evts.t{1}(1,:), 'y'); 
vline(evts.t{1}(2,:), 'y'); 



%% get the movement from the encoder ticks

mov_rate = MS_spike2rate(evts, csc.tvec, 0.1)

%% count the SWR per condition

swr_red = []; 
for ii = length(evts.t{2}):-1:1

 this_data= restrict(swr, evts.t{2}(1,ii), evts.t{2}(2,ii)); 
swr_red(ii) = length(this_data.tstart); 

 this_data= restrict(swr, evts.t{2}(1,ii)-2.5, evts.t{2}(1,ii)); 
ctrl_red(ii) = length(this_data.tstart); 

end

swr_blue = []; 
for ii = length(evts.t{3}):-1:1

 this_data= restrict(swr, evts.t{3}(1,ii), evts.t{3}(2,ii)); 
swr_blue(ii) = length(this_data.tstart); 

end
