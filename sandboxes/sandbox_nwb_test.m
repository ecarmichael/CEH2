addpath('C:\Users\ecar\Documents\Github\matnwb')

cd('C:\Users\ecar\Documents\Github\Toothy\Toothy-main\example_data')

nwb = nwbRead('allenvcn_session829720705_probeF_trimmed.nwb')

%%
electricalseries_map = nwb.searchFor('ElectricalSeries')

all_electricalseries_paths = electricalseries_map.keys();      % Cell array of paths
first_path = all_electricalseries_paths{1};

% Retrieve the object using its path
electricalseries_obj = nwb.resolve(first_path);

% Most data objects have a .data property
tvec = electricalseries_obj.timestamps.load();
raw_data = electricalseries_obj.data.load();
elec = electricalseries_obj.electrodes.data;


% raw_data_size = size(raw_data)

% Check for additional metadata
fprintf('Description: %s\n', electricalseries_obj.description);
%% check the cscs;
figure(101)
clf; 
hold on
for ii  = 1:size(raw_data,1)

    plot(tvec, raw_data(ii,:)*1000 + ii)
end


%% convert to csc data

csc = tsd(tvec, raw_data); 

for ii = 1:size(csc.data,1)

    csc.label{ii} = num2str(ii);
    csc.cfg.hdr{ii}.SamplingFrequency = 1/mode(diff(tvec)); 
end