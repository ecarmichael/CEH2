function c_ord = MS_linspecer(n_colours)
%% MS_linspecer: wrapper for linspecer to remove yellows. 
%
%
%
%    Inputs: 
%    - n_colours: [double]    desired number of colours. 
%
%
%
%    Outputs: 
%    - c_ord:[N x 3 double]   RGB values for colour range. 
%
%
%
%
% EC 2023-11-01   initial version 
%
%
%
%% initialize
y_range = floor(n_colours/10); 
c_ord = linspecer(n_colours+y_range);


c_ord(floor(n_colours/2)-ceil(y_range/2)+1:floor(n_colours/2)+ceil(y_range/2),:) = [];

if size(c_ord,1) < n_colours
    c_ord = linspecer(n_colours+y_range+1);


c_ord(floor(n_colours/2)-ceil(y_range/2)+1:floor(n_colours/2)+ceil(y_range/2),:) = [];
    
    
    
end