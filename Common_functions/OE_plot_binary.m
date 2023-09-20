function OE_plot_binary(stream, channels, range)
%% OE_plot_binary:
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
% EC 2023-09-18   initial version 
%
%
%
%% initialize

if nargin <1
    channels = 1:stream.metadata.numChannels;
    range = [tream.timestamps(1)    stream.timestamps(end)];
elseif nargin <2
        range = [tream.timestamps(1)    stream.timestamps(end)];
end

tstart = nearest_idx3( range(1),stream.timestamps); 
tend = nearest_idx3( range(2),stream.timestamps); 

%% plot

figure
hold on
offset = 1000; 
c_ord = winter(ceil(length(channels)));
c_ord_r = flipud(c_ord); 
for ii = 65:length(channels)
    if rem(ii, 2)
    plot(stream.timestamps(tstart:tend), stream.samples(channels(ii),tstart:tend)+ii*offset, 'color', (c_ord(ii,:)))
    else

    plot(stream.timestamps(tstart:tend), stream.samples(channels(ii),tstart:tend)+ii*offset, 'color',(c_ord_r(ii,:)))
    end


end
