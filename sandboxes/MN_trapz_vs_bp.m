     
csc = MS_LoadCSC([])

    fs = csc.cfg.hdr{1}.SamplingFrequency; % get the sampling frequency from the header. 

win_s = 1; % window size in seconds. 
hann_win = 2^11; %hanning window size
[pxx,f]  = pwelch(csc.data, hanning(hann_win), hann_win/2, hann_win*2 ,fs);
% same thing but with the 'white filtered'
[pxx_d,f]  = pwelch(diff(csc.data), hanning(hann_win), hann_win/2, hann_win*2 ,fs);


%% plot to check

figure(1)
clf
subplot(2,2,1)


plot(f, 10*log10(pxx))
xlim([0 120])

hold on
vline([6 12], {'g' 'g'});
vline([30 58], {'r' 'r'});
vline([62 90], {'b' 'b'});


subplot(2,2,2)

plot(f, 10*log10(pxx_d))
xlim([0 120])

hold on
vline([6 12], {'g' 'g'});
vline([30 58], {'r' 'r'});
vline([62 90], {'b' 'b'});
%% get the AUC using trapizoid
% baseline
base_f = [0 58];
base_f_idx = nearest_idx3( base_f, f);

% theta
theta_f = [6 12];
theta_f_idx = nearest_idx3( theta_f, f);

% low gamma
lg_f = [30 58];
lg_f_idx = nearest_idx3( lg_f, f);

% high gamma
hg_f = [62 90];
hg_f_idx = nearest_idx3(hg_f, f);

% get the AUC for each
base_area = trapz(pxx_d(base_f_idx(1):base_f_idx(2)));
theta_area = trapz(pxx_d(theta_f_idx(1):theta_f_idx(2)));
lg_area = trapz(pxx_d(lg_f_idx(1):lg_f_idx(2)));
hg_area = trapz(pxx_d(hg_f_idx(1):hg_f_idx(2)));

% normalize to the baseline
theta_norm = theta_area/base_area;
lg_norm = lg_area/base_area;
hg_norm = hg_area/base_area;



%% add a bar plot
subplot(2,2,3)

bar(1:3, [theta_norm, lg_norm, hg_norm]')
ylabel('normalized PSD AUC');
set(gca, 'XTickLabel', {'theta', 'lg', 'hg'})



%% same thing but with bandpower

base_bp =  bandpower(pxx_d, f, [0 58], 'psd'); 
theta_bp = bandpower(pxx_d, f, [6 12], 'psd')./ base_bp; 
lg_bp = bandpower(pxx_d, f, [30 58], 'psd')./ base_bp; 
hg_bp = bandpower(pxx_d, f, [62 90], 'psd')./ base_bp; 


%% make a plot of the bandpower

subplot(2,2,4)
b = bar(1:3, [theta_bp, lg_bp, hg_bp]', 'facecolor', 'r')
ylabel('normalized bandpower');
set(gca, 'XTickLabel', {'theta', 'lg', 'hg'})


