function waves = HF_phy2_waveforms(phy_dir, S)
%% HF_phy2_waveforms:
%
%
%
%    Inputs:
%    -
%
%
%
%    Outputs:
%    -
%
%
%
%
% EC 2026-10-28   initial version
%
%
%
%% initialize

if nargin < 2
    error('HF_phy2_waveforms requires phy_dir and S as inputs.');
end
if ~isfolder(phy_dir)
    error('The specified phy_dir does not exist: %s', phy_dir);
end

% parameters

params = OE_load_params(phy_dir);
Fs = params.sample_rate;
nChan = params.n_channels_dat;

nWaves = 5000;
win_s = -round(1/1000*Fs):round(2/1000*Fs);
win_tvec = (win_s/Fs)*1000; % convert time back to miliseconds.
win_samples = numel(win_s);




%% memory mapt the raw data.
% get the location of the raw binary file

params = OE_load_params(phy_dir);

d_idx = strfind(params.dat_path, 'Record');
dat_dir = strrep(strrep([phy_dir(1:strfind(phy_dir, 'kilosort4')-1) params.dat_path(d_idx:end)], '\', filesep), '/', filesep);

ts = readNPY([phy_dir filesep 'timestamps.npy']);
nSamples = size(ts,1);

m = memmapfile(dat_dir, 'Format', {'int16', [nChan, nSamples], 'x'});



%% loop over cells and extract some waveforms.
waves = NaN(length(S.t), win_samples);
w_sign = zeros(length(S.t),1);

for ii = 1:length(S.t)

    if length(S.t{ii}) < nWaves
        this_idx = 1:length(S.t{ii});
    else
        this_idx = datasample(1:length(S.t{ii}), nWaves);
    end

    these_ts = sort(S.t{ii}(this_idx));
    start_idx = nearest_idx3(these_ts, ts-ts(1)) ;
    start_idx(start_idx<1) = [];
    start_idx(start_idx+win_s(end)>length(ts)) = [];

    this_wave = [];
    for jj = length(start_idx) :-1:1
        % for cc = nChan:-1:1
        this_wave(jj,:) = m.Data.x(S.usr.ch(ii)+1, start_idx(jj) + win_s(1):start_idx(jj)+win_s(end));

        this_base = (this_wave(jj,1:abs(win_s(1))));
        this_wave(jj,:) = (this_wave(jj,:) - mean(this_base))/std(this_base);     %(this_wave(jj,:) - this_base)-1;
        % end
    end

    waves(ii,:) = median(this_wave,1);


    % Flip waveforms so their dominant deflection is negative.
    z_idx = nearest_idx(0, win_tvec); 


    %check if the trough is before the zero point and if so, redo the  
   [all_t all_t_idx] = min(waves(ii,:)); 

   % if all_t_idx < z_idx
   % 
   %     this_wave = [];
   %     for jj = length(start_idx) :-1:1
   %         % for cc = nChan:-1:1
   %         this_wave(jj,:) = m.Data.x(S.usr.ch(ii)+1, start_idx(jj) + win_s(1) - (z_idx-all_t_idx):start_idx(jj)+win_s(end)-(z_idx-all_t_idx));
   % 
   %         this_base = (this_wave(jj,1:abs(win_s(1))));
   %         this_wave(jj,:) = (this_wave(jj,:) - mean(this_base))/std(this_base);     %(this_wave(jj,:) - this_base)-1;
   %         % end
   %     end
   % 
   %     waves(ii,:) = median(this_wave,1);
   % 
   % end



    [t, t_idx] = min((waves(ii,z_idx:end)));
    [p, p_idx] = max((waves(ii,1:(t_idx+z_idx)-1)));
    [p_all, p_all_idx] = max(waves(ii,:)); 

    if abs(p)>abs(t)  || t > 0  % p_all_idx < all_t_idx
    %if p_idx < t_idx+z_idx
        % if waves(ii, nearest_idx(0, win_tvec)) > 0
        waves(ii,:) = -waves(ii,:);
        w_sign(ii) = 1; 
    end

    all_waves{ii} = this_wave;

end

%%
S_out = MS_get_wave_properties_phy2(S, waves,win_tvec,  1);



%% test umap

c_ord = MS_linspecer(3); 

wave_prop = [S_out.usr.dur; S_out.usr.pt_ratio; S_out.usr.pre_t_ratio;S_out.usr.pre_t_dur; S_out.usr.pre_peak_ratio; S_out.usr.repol_slope; S_out.usr.recov_slope]';

man_class = [0 1 0 1 0 1 1,...
             1 1 1 1 1 1 1,...
             0 0 0 0 1 1 1,...
             1 1 0 1 0 0 0,...
             0 0 0 0 0 0 0,...
             0 0 0];

wave_map = umap(wave_prop, 'NumDimensions',3);

figure(900)
clf
subplot(2,2,1)

scatter(wave_prop(:,1),wave_prop(:,2), 50,c_ord(man_class+1,:), 'filled')
xlabel('duration'); 
ylabel('pt ratio')
zlabel('pre t ratio')

subplot(2,2,2)
scatter(wave_map(:,1),wave_map(:,2), 50,c_ord(man_class+1,:), 'filled')

title("2-D Embedding")%% wave characteristics based on NeuroPix Ultra paper

subplot(2,2,3)

scatter3(wave_prop(:,1),wave_prop(:,2), wave_prop(:,3), 50,c_ord(man_class+1,:), 'filled')
xlabel('duration'); 
ylabel('pt ratio')
zlabel('pre t ratio')

subplot(2,2,4)
scatter3(wave_map(:,1),wave_map(:,2), wave_map(:,3), 50,c_ord(man_class+1,:), 'filled')

title("3-D Embedding")%% wave characteristics based on NeuroPix Ultra paper




%% quick plot
figure(101)
clf
pause(.1)
n = ceil(length(S.t)/7);
mm = ceil(length(S.t)/6);

for cc = 1:length(S.t)
    subplot(n,mm,cc)
    hold on
        plot(win_s,waves(cc,:), 'k')

    if w_sign(cc)
        plot(win_s, -waves(cc,:), '--r');
    else
        plot(win_s, waves(cc,:),  '--r');
    end
    xline(0, 'k:');
    title(sprintf('Unit %d', cc));

end



end


% %% read the phy spikes waves
%
% wave_chan = readNPY('_phy_spikes_subset.channels.npy');
% wave_data = readNPY('_phy_spikes_subset.waveforms.npy');
% wave_spikes = readNPY('_phy_spikes_subset.spikes.npy');
%
% wave_data = reshape(wave_data, [], win_samples, size(wave_data, 3));
% wave_chan = wave_chan(:);
% waves = squeeze(median(wave_data, 1, 'omitnan'));