function [stream] = OE_load_binary(data_dir); 
%% OE_load_binary: quickway to load binary format OE data. Just pieces I grabed from BinaryRecording in https://www.mathworks.com/matlabcentral/fileexchange/122372-open-ephys-matlab-tools
%
%
%
%    Inputs: 
%    - data_dir[path] path to continuous.dat file and
%    timestamps/sample_numbers. 
%
%
%
%    Outputs: 
%    - data_out:[struct]  data in the TSD format. 
%
%
%
%
% EC 2023-09-18   initial version 
%
%
%
%% initialize

cd(data_dir)

up_dir = fullfile(cd, '..');
cd(up_dir)
up_dir = fullfile(cd, '..');
cd(up_dir)

self.info = jsondecode(fileread(fullfile(cd,'structure.oebin'))); 
self.info.continuous(1);

 stream.metadata.sampleRate = self.info.continuous(1).sample_rate;
                stream.metadata.numChannels = self.info.continuous(1).num_channels;
                stream.metadata.processorId = self.info.continuous(1).source_processor_id;
                stream.metadata.streamName = self.info.continuous(1).folder_name(1:end-1);

%%
cd(data_dir)

stream.timestamps = readNPY(fullfile(data_dir, 'timestamps.npy'));

data = memmapfile(fullfile(data_dir, 'continuous.dat'), 'Format', 'int16');
stream.metadata.numChannels = 67;
stream.samples = reshape(data.Data, [stream.metadata.numChannels, length(data.Data) / stream.metadata.numChannels]);
stream.sampleNumbers = readNPY(fullfile(data_dir, 'sample_numbers.npy'));

%% refomat



