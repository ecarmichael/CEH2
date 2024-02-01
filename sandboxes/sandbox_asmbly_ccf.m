figure(1010)
clf

w_data = wake_xcor{iB}; 
Pre_data = REM_pre_xcor{iB};
Post_data = REM_post_xcor{iB}; 

w_dataz = wake_zxcor{iB}; 
Pre_dataz = REM_pre_zxcor{iB};
Post_dataz = REM_post_zxcor{iB}; 

sig_idx = wake_zxcor{iB} > 1.96; 

Pre_sig_mat = Pre_dataz > 1.96;


c_lim = [1.96 inf]; 

% raw cor
subplot(2,3,1)
cla
imagesc(Pre_data)
title('pre xcorr')
% caxis(c_lim);

subplot(2,3,2)
cla
imagesc(w_data)
title('wake xcorr')

% caxis(c_lim);


subplot(2,3,3)
cla
imagesc(Post_data)
title('post xcorr')

% caxis(c_lim);

%
subplot(2,3,4)
cla
imagesc(Pre_dataz)
title('pre xcorr Z')
caxis(c_lim);

subplot(2,3,5)
cla
imagesc(w_dataz)
caxis(c_lim);
title('wake xcorr Z')


subplot(2,3,6)
cla
imagesc(Post_dataz)
caxis(c_lim);
title('post xcorr Z')


figure(1011)

subplot(1,2,1)
cla
imagesc(Pre_data./ w_data)
title('pre/wake')

subplot(1,2,2)
cla
imagesc(Post_data ./w_data)
title('post/wake')

% print the number of signficant ccf pairs from the awake sig pairs

fprintf('Pre REM had %.0f (%.0f%%) significant assemblies of the %0.0f found in wake \n',...
    (sum(Pre_dataz > 1.96, 'all') - numel(diag(sig_idx))), ((sum(Pre_dataz > 1.96, 'all') - numel(diag(sig_idx))) / (numel(sig_idx) - numel(diag(sig_idx))))*100,...
    (sum(sig_idx, 'all') - numel(diag(sig_idx))))

fprintf('Post REM had %.0f (%.0f%%) significant assemblies of the %0.0f found in wake \n',...
    (sum(Post_dataz > 1.96, 'all') - numel(diag(sig_idx))), ((sum(Post_dataz > 1.96, 'all') - numel(diag(sig_idx))) / (numel(sig_idx) - numel(diag(sig_idx))))*100,...
    (sum(sig_idx, 'all') - numel(diag(sig_idx))))


%% generate fake assemblies.

shift_data = [];
for ii =1:15
        
    shift_data(:,ii) = circshift(data_h(:,1), ii);
    
    
    
end

data_toy = [data_toy shift_data]; 



%% csp test plots
figure(10102)
iC = 3; 
clf
hold on
        opts.binsize = 1/30; 
        opts.gauss_window = 1./opts.binsize; % 1 second window
        opts.gauss_SD = 0.1./opts.binsize; % 0.02 seconds (20ms) SD
        gk = gausskernel(opts.gauss_window,opts.gauss_SD); gk = gk./opts.binsize; % normalize by binsize

        

    Csp =all_deconv_pre_REM./all_denoise_pre_REM;
    Csp(isnan(Csp)) = 0; 
    Csp = Csp > 0.01;
    REM_pre_data = Csp;


t = (0:size(REM_pre_data,1)-1)./30; 
plot(t, REM_pre_data(:,iC), 'r')

plot(t, all_deconv_pre_REM(:,iC)./all_denoise_pre_REM(:,iC), 'm')
plot(t, all_deconv_pre_REM(:,iC)*5, 'k')

plot(t, conv2(REM_pre_data(:,iC), gk, 'same'), 'g')
% 


%% test the MAD method
iC =500
figure(909)
clf
trace = ms.detrendRaw(:,iC);
d_trace = ms.denoise(:,iC);

ax(1) = subplot(4,1,1:3);
hold on
plot(ms.time, trace, 'r')
plot(ms.time, d_trace, 'k')

ax(2) = subplot(4,1,4);
hold on
plot(ms.time, ms.deconv(:,iC)./m_ad, 'b')

m_ad = mad([d_trace- trace]);

thresh = 1.5; 
yline(thresh)
linkprop(ax, 'xlim')
keep_idx = ms.deconv(:,iC) > thresh; 
plot(ms.time(keep_idx), ms.deconv(keep_idx,iC)./m_ad, '*m')


