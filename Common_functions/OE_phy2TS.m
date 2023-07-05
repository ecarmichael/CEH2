function S = OE_phy2TS()

% tvec = readNPY('timestamps.npy');

spike_times = readNPY('spike_times.npy');
spike_clusters = readNPY('spike_clusters.npy');

fid = fopen('cluster_group.tsv');
C = textscan(fid, '%s%s');
fclose(fid);

fid = fopen('cluster_info.tsv');
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
    S.t{ii} = double(spike_times(spike_clusters == good_clusters_ids(ii)));
    S.label{ii} = [num2str(good_clusters_ids(ii)) '-' num2str(I(ii))]; 

end

S.cfg.history.mfun{1} = mfilename;
S.cfg.history.cfg{1} = [];

end