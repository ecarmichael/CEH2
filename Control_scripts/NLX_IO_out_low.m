function [succeeded, CheetahReply] = NLX_IO_out_low(port, pin)
% Initialize the output variables
succeeded = false;
CheetahReply = [];

% sent the command
[succeeded,CheetahReply]=NlxSendCommand(cat(2,'-DigitalIOTtlPulse AcqSystem1_0 ', num2str(port), ' ', num2str(pin), ' Low'));
