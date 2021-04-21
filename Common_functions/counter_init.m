function counter_init(iter, outof);
%% counter_init is this initiation counter.   
%
%
%
%   Inputs: 
%       - iter: value to print. Normally your iteration. [required]
%
%       - outof: denominator. [optional]
%
%
%
%
%
%
% EC 2021-04-19   initial version 
%
%
%
%% initialize
if nargin ==1   
    fprintf(1,'\n%0.0f...',iter) % Deleting 4 characters (The three digits and the % symbol)
else
    fprintf(1,'\n%0.0f/%0.0f...',iter, outof) % Deleting 4 characters (The three digits and the % symbol)
end