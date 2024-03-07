function OE_theta_test(data_dir, chan, evt_id)


%% find the data and load it. 

if isempty(data_dir)
    data_dir = uigetdir;
end

%% load the csc
csc_dir = dir([data_dir '\**\*.continuous']);

for ii = 1:length(csc_dir)
   
    if contains(csc_dir(ii).name, ['CH' chan])
        keep_idx(ii) = 1; 
    else
        keep_idx(ii) = 0; 
    end
end

fprintf('Using channel <strong>%s</strong>...\n', csc_dir(find(keep_idx)).name)


cfg.fc = {csc_dir(find(keep_idx)).name}; 
cfg.desired_sampling_frequency = 2000;
csc = OE_old_csc2TSD(cfg);

% load some events
evts_dir = dir([csc_dir(find(keep_idx)).folder '\*Rhythm*.events']);
evts = OE_LoadEvents(evts_dir.name);

evt_idx = find(contains(evts.label, num2str(evt_id))); 

%% split the TTL times into start and stop


t_start = evts.t{evt_idx}(1:2:end);
t_end = evts.t{evt_idx}(2:2:end);

ttl_iv = iv(t_start, t_end); 

%% compute the psd for each on periods (assuming you started with on
hann_win = 2^12; 

for ii = 1:2:length(evts.t{evt_idx})
    
   csc_r = restrict(csc, evts.t{evt_idx}(ii), evts.t{evt_idx}(ii+1)); 
   
    [Pxx_on{ii}, f] = pwelch(csc_r.data(1,:), hanning(hann_win), hann_win/2, hann_win*2 , csc_r.cfg.hdr{1}.SamplingFrequency);

end

Pxx_on(cellfun(@isempty, Pxx_on)) = []; 
Pxx_on = cell2mat(Pxx_on); 

% same for off
for ii = 2:2:length(evts.t{evt_idx})-1
    
   csc_r = restrict(csc, evts.t{evt_idx}(ii), evts.t{evt_idx}(ii+1)); 
   
    [Pxx_off{ii}, f] = pwelch(csc_r.data(1,:), hanning(hann_win), hann_win/2, hann_win*2 , csc_r.cfg.hdr{1}.SamplingFrequency);

end

Pxx_off(cellfun(@isempty, Pxx_off)) = []; 
Pxx_off = cell2mat(Pxx_off); 


%% get the mean PSD

figure(101)
clf
subplot(4,2,1:4)
plot(csc)
hold on
for ii = 1:2:length(evts.t{evt_idx})
   rectangle('Position', [ evts.t{evt_idx}(ii), max(csc.data), evts.t{evt_idx}(ii+1) - evts.t{evt_idx}(ii), max(csc.data)/20], 'FaceColor', 'g', 'EdgeColor', 'g')
   
%    evts.t{evt_idx}(ii):evts.t{evt_idx}(ii+1), max(csc.data)*ones(length([evts.t{evt_idx}(ii):evts.t{evt_idx}(ii+1)]),1), 'g')
    
end

subplot(4,2,5)

hold on
plot(f, Pxx_on', 'g', 'linewidth', 1)

xlim([0 20]);ylim([0 max([Pxx_on Pxx_off],[], 'all')]);
xlabel('Freq (Hz)')
ylabel('Power (db)')

subplot(4,2,6)

hold on
plot(f, Pxx_off', 'color', [.7 .7 .7], 'linewidth', 1)

xlim([0 20]);ylim([0 max([Pxx_on Pxx_off],[], 'all')]);xlabel('Freq (Hz)')
ylabel('Power (db)')





subplot(4,2,7:8)

hold on
plot(f, mean(Pxx_on'), 'g', 'linewidth', 3)
plot(f, mean(Pxx_off'), 'color', [.7 .7 .7], 'linewidth', 3)

xlim([0 20])
xlabel('Freq (Hz)')
ylabel('Power (db)')

