%% RunClustBatch_M26

% TT1

if exist(['FD' filesep 'TT1.clu.1']); delete(['FD' filesep 'TT1.clu.1']);end
RunClustBatch('fcTT', {'TT1.ntt'}, 'channelValidity',  [1 1 1 1], 'minClusters', 10, 'maxClusters', 20,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

% TT2
if exist(['FD' filesep 'TT2.clu.1']); delete(['FD' filesep 'TT2.clu.1']);end
RunClustBatch('fcTT', {'TT2.ntt'}, 'channelValidity',  [1 0 1 1], 'minClusters', 10, 'maxClusters', 20,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

% TT3
if exist(['FD' filesep 'TT3.clu.1']); delete(['FD' filesep 'TT3.clu.1']);end
RunClustBatch('fcTT', {'TT3.ntt'}, 'channelValidity',  [1 0 0 1], 'minClusters', 10, 'maxClusters', 20,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});
% 
% TT4
if exist(['FD' filesep 'TT4.clu.1']); delete(['FD' filesep 'TT4.clu.1']);end
RunClustBatch('fcTT', {'TT4.ntt'}, 'channelValidity',  [1 1 1 1], 'minClusters', 10, 'maxClusters', 20,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

% TT5

if exist(['FD' filesep 'TT5.clu.1']); delete(['FD' filesep 'TT5.clu.1']);end
RunClustBatch('fcTT', {'TT5.ntt'}, 'channelValidity',  [0 1 0 1], 'minClusters', 10, 'maxClusters', 20,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

% TT6
if exist(['FD' filesep 'TT6.clu.1']); delete(['FD' filesep 'TT6.clu.1']);end
RunClustBatch('fcTT', {'TT6.ntt'}, 'channelValidity',  [1 1 1 1], 'minClusters', 10, 'maxClusters', 20,...
    'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});
% 
% % TT7
% 
% if exist(['FD' filesep 'TT7.clu.1']); delete(['FD' filesep 'TT7.clu.1']);end
% RunClustBatch('fcTT', {'TT7.ntt'}, 'channelValidity',  [0 1 1 1], 'minClusters', 10, 'maxClusters', 20,...
%     'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});
% 
% % TT8
% if exist(['FD' filesep 'TT8.clu.1']); delete(['FD' filesep 'TT8.clu.1']);end
% RunClustBatch('fcTT', {'TT8.ntt'}, 'channelValidity',  [1 1 1 0], 'minClusters', 10, 'maxClusters', 20,...
%     'featuresToCompute',{'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'}, 'featuresToUse', {'feature_Peak','feature_Valley','feature_Energy', 'feature_EnergyD1'});

