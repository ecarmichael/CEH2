%% RunClustBatch_M26

% TT1
RunClustBatch('fcTT', {'TT1.ntt'}, 'channelValidity',  [0 0 1 1], 'minClusters', 6, 'maxClusters', 12,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

% TT2
RunClustBatch('fcTT', {'TT2.ntt'}, 'channelValidity',  [1 0 1 1], 'minClusters', 6, 'maxClusters', 12,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

% TT3
RunClustBatch('fcTT', {'TT3.ntt'}, 'channelValidity',  [1 0 0 1], 'minClusters', 6, 'maxClusters', 12,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

% TT4
RunClustBatch('fcTT', {'TT4.ntt'}, 'channelValidity',  [1 1 1 1], 'minClusters', 6, 'maxClusters', 12,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

% TT5
RunClustBatch('fcTT', {'TT5.ntt'}, 'channelValidity',  [0 1 0 1], 'minClusters', 6, 'maxClusters', 12,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

% TT6
RunClustBatch('fcTT', {'TT6.ntt'}, 'channelValidity',  [1 1 1 1], 'minClusters', 6, 'maxClusters', 12,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

% TT7
RunClustBatch('fcTT', {'TT7.ntt'}, 'channelValidity',  [1 1 1 1], 'minClusters', 6, 'maxClusters', 12,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

% TT8
RunClustBatch('fcTT', {'TT8.ntt'}, 'channelValidity',  [1 1 1 0], 'minClusters', 6, 'maxClusters', 12,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

