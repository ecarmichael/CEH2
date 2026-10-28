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

nWaves = 10000; 
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
    if waves(ii,abs(win_s(1))) >0
        waves(ii,:) = -waves(ii,:);
    end

    all_waves{ii} = this_wave;

end

%%
for ii = 1:length(S.t)
S_t = SelectTS([], S, ii);
MS_get_wave_properties_phy2(S_t, waves(ii,:),win_tvec,  1)
end


%% wave characteristics based on NeuroPix Ultra paper

wave = waves(3,:);

%zero
z_idx = nearest_idx(0, win_tvec); 

% trough
[t, t_idx] = min(wave);

% (post-trough) peak
[p, p_idx] = max(wave(z_idx:end)); 

% pre peak
[pre, pre_idx] = max(wave(1:z_idx)); 

% duration
dur = win_tvec(p_idx) - win_tvec(t_idx); 

% amplitude

% trough - peak ratio
troughPeakRatio = abs(t / p);

% repolarization slope
repolSlope = (p - t) / (win_tvec(p_idx) - win_tvec(t_idx));

% pre-peak to trough ratio
prePeakTroughRatio = abs(pre / t);
% trough - peak ratio

% repol slope

% prePeak Trough ratio

% flip if needed
% if min(wave(1:z_idx)) > min(wave(z_idx:end))



    %% quick plot
    figure(101)
    clf
    for cc = 1:25
        subplot(5,5,cc)
           plot(win_s,waves(cc,:))
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