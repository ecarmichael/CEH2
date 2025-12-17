function S = OE_phy2TS(data_dir)
%% OE_phy2TS: loads the classified cells from Phy2 and kilosort. 
%
%
%
%   Inputs:
%      - data_dir: [directory with the kilosort / phy outputs
%
%      - label: [string] or [cell array with strings] cell labels to
%      include ('good' 'mua' 'unsorted'). NOT implemented yet. Only good 
%
%
%    Outputs:
%       - S [struct]  spike data in the TS format. contains spike times. 

if nargin < 1
    data_dir = cd; 
end

% tvec = readNPY('timestamps.npy');
spike_struct = loadParamsPy([data_dir filesep 'params.py']);  % from https://github.com/cortex-lab/spikes/blob/master/preprocessing/phyHelpers/loadKSdir.m

spike_times = readNPY([data_dir filesep 'spike_times.npy']);
% convert to time. 
spike_times = double(spike_times)/spike_struct.sample_rate; % from  https://github.com/cortex-lab/spikes/blob/master/preprocessing/phyHelpers/loadKSdir.m
spike_clusters = readNPY([data_dir filesep 'spike_clusters.npy']);

%chan_shanks = readNPY([data_dir filesep 'channel_shanks.npy']);
%chan_pos= readNPY([data_dir filesep 'spike_positions.npy']);


fid = fopen([data_dir filesep 'cluster_group.tsv']);
C = textscan(fid, '%s%s');
fclose(fid);

fid = fopen([data_dir filesep 'cluster_info.tsv']);
I= textscan(fid, '%s%s%s%s%s%s');
fclose(fid);

I =I{6}(3:2:end); 
I = str2double(I); 


cids = cellfun(@str2num, C{1}(2:end), 'uni', false);
ise = cellfun(@isempty, cids);
cids = [cids{~ise}];

isUns = cellfun(@(x)strcmp(x,'unsorted'),C{2}(2:end));
isMUA = cellfun(@(x)strcmp(x,'mua'),C{2}(2:end));
isGood = cellfun(@(x)strcmp(x,'good'),C{2}(2:end));
cgs = zeros(size(cids));

cgs(isMUA) = 1;
cgs(isGood) = 2;
cgs(isUns) = 3;

%% extract the good units
good_clusters_ids = cids(isGood);
S = [];
S.type = 'ts';
for ii = length(good_clusters_ids):-1:1
    S.t{ii} = (spike_times(spike_clusters == good_clusters_ids(ii)));
    S.label{ii} = [num2str(good_clusters_ids(ii)) '-' num2str(I(ii))]; 

end

S.cfg.history.mfun{1} = mfilename;
S.cfg.history.cfg{1} = [];

end