%% sandbox_CIHR_HR_TFC


kilo_dir = []; 

ctrl_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/CIHR_2025/HF/conditioning 1/arduino_log_1756923023.csv'; 

eye_dir = []; 



%% load the control log

log_tab = readtable(ctrl_dir); 

log_tab.time = log_tab.time - log_tab.time(1);  % zero to the first time point. 

log_tab.time = log_tab.time/ 1000; % convert from ms to seconds;

% convert to evts

labels = unique(log_tab.phase); 
phases = {'baseline', 'Tone', 'Trace', 'Puff'}; 

for ii = 1:length(phases)

log.(phases{ii}) = log_tab.time(contains(log_tab.event, '-') & contains(log_tab.phase, phases{ii})); 

end
log.end = log_tab.time(contains(log_tab.event, 'Finish')); 

% get the encoder changes as a vector. 

wheel.tvec = log_tab.time(contains(log_tab.event, 'pos')); 
wheel.data = log_tab.encoderCount(contains(log_tab.event, 'pos')); 

% interp
tvec_i = 0:.01:wheel.tvec(end); 

wheel.data = interp1(wheel.tvec, wheel.data, tvec_i); 
wheel.tvec = tvec_i; 


%% plot??
c_ord = MS_linspecer(8); 

figure(101)

clf

subplot(4,1,1)
cla

plot(wheel.tvec, wheel.data, 'k');

for ii = 1:length(phases)
xline(log.(phases{ii}), 'Color',c_ord(ii,:))
end
ylabel('Encoder pos')

subplot(4,1,2)
cla
plot(wheel.tvec(1:end-1)  - wheel.tvec(1), abs(diff(wheel.data)), 'r');
ylabel('Movement (a.u.)')

for ii = 1:length(phases)
xline(log.(phases{ii}), 'Color',c_ord(ii,:))
end

