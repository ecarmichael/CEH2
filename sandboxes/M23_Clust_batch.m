%% RunClustBatch_M23

% TT1
RunClustBatch('fcTT', {'TT1.ntt'}, 'channelValidity',  [1 1 1 1], 'minClusters', 10, 'maxClusters', 20,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

% % TT2
% RunClustBatch('fcTT', {'TT2.ntt'}, 'channelValidity',  [0 1 1 1], 'minClusters', 10, 'maxClusters', 20,...
%     'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});
% 
% % TT3
% RunClustBatch('fcTT', {'TT3.ntt'}, 'channelValidity',  [1 0 1 1], 'minClusters', 10, 'maxClusters', 20,...
%     'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

% TT4
RunClustBatch('fcTT', {'TT4.ntt'}, 'channelValidity',  [0 1 1 0], 'minClusters', 10, 'maxClusters', 20,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

