%% sandbox_Buz_lin_summary



data_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/test_data/2025-07-11_15-46-53';

cd(data_dir)


%% load all the spikes from the recording

% S = MS_Load_NTT([], 'TT*_0001.ntt'); 
S = MS_Load_NTT([])
% grab the LFP for the linear shank

cfg_lin = []; 
cfg_lin.desired_sampling_frequency = 2000; 
cfg_lin.fc = []; 
for ii = 26:41
    cfg_lin.fc{end+1} = ['CSC' num2str(ii) '.ncs']; %_0001.ncs']; 
end

csc = MS_LoadCSC(cfg_lin); 
%% found the cells per tetrode
TT_count = zeros(1,12); 
for ii = 1:12
    for iS = 1:length(S.label)
        if contains(S.label{iS}, ['TT' num2str(ii) '_'])
            TT_count(ii) = TT_count(ii)+1; 
        end
    end
    fprintf('<strong>%s</strong>: %.0f cells\n', ['TT' num2str(ii)], TT_count(ii))
end


%% compute the autocorrelation for each
AC = cell(size(S.t));
AC_x  =AC;

figure(2)
clf
for ii = 1:length(S.t)
    for jj = 1:length(S.t)
        [AC{ii, jj}, AC_x{ii, jj}] =xcorr(S.t{ii}, S.t{jj}, 100);
  
        subplot(size(S.t,2), size(S.t,2), (ii+((ii-1)*size(S.t,2))))
        plot(AC_x{ii, jj}, AC{ii, jj})
    end


end


%%  compute the CSD around a small window

cfg_filter = [];
cfg_filter.f = [6 12];
cfg_filter.type = 'cheby1';
cfg_filter.order = 3;
% cfg_filter.display_filter = 1;

csc_f = FilterLFP(cfg_filter, csc);

[csd, csd_inds]= csd5pt_Schomburg(csc_f.data, 100);  


cfg_filter = [];
cfg_filter.f = [30 55];
cfg_filter.type = 'butter';
cfg_filter.order = 3;
% cfg_filter.display_filter = 1;

csc_g = FilterLFP(cfg_filter, csc);
[csd_g, csd_inds]= csd5pt_Schomburg(csc_g.data, 100);  

%%
figure(1)
clf
hold on
% imagesc(csc_f.tvec - csc_f.tvec(1),1:size(csd,1),csd)

for ii = 1:length(csd_inds)
    plot(csc_f.tvec - csc_f.tvec(1), (csc.data(csd_inds(ii),:)*3500)+ii, 'k');
    plot(csc_f.tvec - csc_f.tvec(1), (csc_f.data(csd_inds(ii),:)*1500)+ii, 'b');
    % plot(csc_g.tvec - csc_g.tvec(1), abs(hilbert((csc_g.data(csd_inds(ii),:)))*2000)+ii, 'r');
    plot(csc_g.tvec - csc_g.tvec(1), (csc_g.data(csd_inds(ii),:)*2000)+ii, 'r');


end

set(gca, 'ytick', 1:length(csd_inds), 'YTickLabel', csc.label(csd_inds), 'YDir', 'normal')


%% average over theta cycles


csc_f_r = restrict(csc_f, csc.tvec(1)+205, csc.tvec(1)+226);


[csd_r, csd_inds]= csd5pt_Schomburg(csc_f_r.data, 100);  


phi = angle(hilbert(csc_f_r.data(6,:)));

phi_peaks = []; 
 for ii = 1:length(phi)-1
        if (phi(ii)>0) && (phi(ii+1)<=0)
            phi_peaks(end+1) = ii+1;
        end
 end

 all_t_cycle = []; all_csd_cycle = []; 
 for ii = length(phi_peaks)-10:-1:10
     all_t_cycle(ii,:) = csc_f_r.data(1,phi_peaks(ii)-500:phi_peaks(ii)+500);
     all_csd_cycle(:,:,ii) = csd_r(:,phi_peaks(ii)-500:phi_peaks(ii)+500); 
 end

  all_t_cycle = all_t_cycle(10:length(phi_peaks)-10,:); 

 all_csd_cycle = all_t_cycle(:,:,10:length(phi_peaks)-10); 

 figure(1818)
 clf
 hold on
 plot(csc_r.tvec - csc_r.tvec(1), csc_f_r.data(1,:), 'b')
 plot(csc_r.tvec - csc_r.tvec(1), csc_r.data(1,:), 'k')
 plot(csc_r.tvec - csc_r.tvec(1), phi/10000, 'r')
 plot(csc_r.tvec(phi_peaks) - csc_r.tvec(1), (pi/2)/10000, 'xm')

%% Phase amp mod
csc_r = restrict(csc, csc.tvec(1)+205, csc.tvec(1)+226);

c = 6; 
c_idx = zeros(size(csc_r.label)); 
c_idx(6) = 1; 

csc_p = csc_r; 
csc_p.data(~c_idx,:) = [];  
csc_p.label(~c_idx) = []; 

cfg_m.P_step = .2; % step size in Hz
cfg_m.A_step = 1;
[como, p, a] = MS_phase_freq(cfg_m, csc, [4 12], [30 120] ); 

figure(101)
clf
hold on
imagesc(p, a, como')
plot(mean(como,2))



%% 

c = 1; 
c_idx = zeros(size(csc.label)); 
c_idx(c) = 1; 

csc_q = csc; 
csc_q.data(~c_idx,:) = [];  
csc_q.label(~c_idx) = []; 

csc_q.tvec = csc_q.tvec-csc_q.tvec(1); 


S_q = S;
for ii = 1:length(S.t)
    S_q.t{ii} = S.t{ii} -csc.tvec(1); 
end

cfg_plt = []; 
cfg_plt.lfp = csc_q; 

MultiRaster(cfg_plt, S_q)
ylim([-20, length(S.label)])

%% 

swr_evts = MS_SWR_detector(csc)