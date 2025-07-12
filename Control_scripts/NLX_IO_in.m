function [trigger] = NLX_IO_in(port)
% Initialize the output variables
trigger = [];

% sent the command
[~,CheetahReply]=NlxSendCommand(cat(2,'-GetDigitalIOPortString AcqSystem1_0 ', num2str(port)));

trigger = bin2dec(CheetahReply);