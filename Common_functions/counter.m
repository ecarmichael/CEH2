function counter(iter, outof)
%% counter just fprintf's a counter that deletes after each iteration.  
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
% if nargin ==1   
%     fprintf([repmat('\b', 1, strlength([num2str(iter) '...'])) '%0.0f...'],iter) % Deleting 4 characters (The three digits and the % symbol)
% else
%     fprintf([repmat('\b', 1, strlength([num2str(iter) '/' num2str(outof) '...'])) '%0.0f/%0.0f...'],iter, outof) % Deleting 4 characters (The three digits and the % symbol)
% end


	fprintf(1,'\b\b\b\b%3.0f%%',100*(iter/outof)); 