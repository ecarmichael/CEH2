function vals =  MS_randn_range(n, m, min_val, max_val)
%% MS_randn_range: get a set of random numbers n x m matrix ranging between min_val and max_val. 
%
%   Example:   vals = MS_randn_range(1,100,0,5)
%
%    Inputs: 
%    - n: number Rows
%
%    - m: number of columns
%
%    - min_val: lowest value for range
%
%    - max_val: largest value
%    Outputs: 
%    - vals [n x m]  random values between min_val and max_val. 
%
%
% EC 2020-04-30   initial version 
%
%
%
%% get the values

vals = (max_val-min_val).*rand(n,m) + min_val;

%fprintf('<strong>%s</strong>: min: %d  max: %d\n', mfilename, min(vals), max(vals));