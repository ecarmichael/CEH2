function tsd_out = MS_behav2tsd(behav, conv_fact)
%% MS_behav2tsd: convert behav to tsd format. 
%
%
%
%   Input:
%       behav:  [struct] contains time and position data
%
%
%
%   output:
%       tsd_out: [struct] contains position 

%% convert everything
if nargin<2
    conv_fact = 1;
    tsd_out = tsd(behav.time, behav.position(:,1:2)'/conv_fact, {'x', 'y'}, 'px');
else
    tsd_out = tsd(behav.time, behav.position(:,1:2)'/conv_fact, {'x', 'y'}, 'cm');
end
