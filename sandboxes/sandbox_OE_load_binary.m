function OE_load_binary(data_dir)
%% OE_load_binary:
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
% EC 2023-08-31   initial version 
%
%
%
%% initialize

session = Session(data_dir);
node = session.recordNodes{1};
recording= node.recordings{1,1};
streamnames = recording.continuous.keys();
streamName = streamnames{1}
data = recording.continuous(streamName);

%% plot the continuous data

lim = 30000*2; % in samples

figure(1)
clf
hold on

offset = 10; 
ticks =[];
for ii = 1:64
    if connected(ii)
    
    plot(data.timestamps(1:lim), (data.samples(ii,1:lim)/100)+(ii*offset))
    ticks(ii) = median((data.samples(ii,1:3000)/100)+(ii*offset)); 
    
    else
            plot(data.timestamps(1:lim), (data.samples(ii,1:lim)/100)+(ii*offset), 'color', [0.8 0.8 0.8])
    ticks(ii) = median((data.samples(ii,1:3000)/100)+(ii*offset)); 
        
    end
end


set(gca, 'YTick', ticks, 'YTickLabel', num2str(chanMap))

% ylim([ticks(1)-ticks(1)*0.1, ticks(end) + ticks(end)*10])

