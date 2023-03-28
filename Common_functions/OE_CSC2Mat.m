function OE_CSC2Mat(fname)
%% OE_CSC2Mat:
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
% EC 2023-03-24   initial version 
%
%
%
%% initialize


session = Session(cd);
node = session.recordNodes{1};

rec = node.recordings{1,1}; 
streams = rec.continuous.keys();

data = rec.continuous(streams{1});
csc = data.samples; 
tvec = data.timestamps; 
save('CSC', 'csc');
save('tvec', 'tvec');