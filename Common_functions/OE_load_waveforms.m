function waves = OE_load_waveforms(phy_dir)
error('not working ATM')


disp('here')

spike_struct = loadParamsPy([phy_dir filesep 'params.py']);  % from https://github.com/cortex-lab/spikes/blob/master/preprocessing/phyHelpers/loadKSdir.m

spike_ind = readNPY([phy_dir filesep 'spike_times.npy']);

% convert to time. 
spike_times = double(spike_ind)/spike_struct.sample_rate; % from  https://github.com/cortex-lab/spikes/blob/master/preprocessing/phyHelpers/loadKSdir.m
spike_clusters = readNPY([phy_dir filesep '_phy_spikes_subset.channels.npy']);
spike = readNPY([phy_dir filesep '_phy_spikes_subset.spikes.npy']);
spike_waves = readNPY([phy_dir filesep '_phy_spikes_subset.waveforms.npy']);


spks = unique(spike_clusters(:,1));
waves = NaN(length(spks), size(spike_waves,2), 4); 

for ii = 1:length(spks)
    
    this_idx = spike_clusters(:,1) == spks(ii); 

    this_wave = spike_waves(this_idx, :, 1:4);
    waves(:,:, 1:4) = spike_waves(this_idx, :, 1:4);
    

end