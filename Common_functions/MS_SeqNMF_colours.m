function S_col = MS_SeqNMF_colours(K)
%% MS_SeqNMF_colours: just a way to pay homage to the SeqNMF colour scheme. from https://github.com/FeeLab/seqNMF/blob/master/SimpleWHPlot.m   in https://elifesciences.org/articles/38471
%
%
%
%    Inputs: 
%    - K: number of colors to pick. limited to 15 if you want to use the OG
%    scheme. 
%
%
%
%    Outputs: 
%    - S_col [K x 3] colors in RGB. 
%
%
%
%
% EC 2022-04-25   initial version 
%
%
%
%% initialize

S_col = [[0 .6 .3]; [.7 0 .7]; [1 .6 0];  [.1 .3 .9];  [1 .1 .1];  [0 .9 .3]; [.4 .2 .7]; [.7 .2 .1]; [.1 .8 1 ]; [1 .3 .7]; [.2 .8 .2]; [.7 .4 1]; [.9 .6 .4]; [0 .6 1]; [1 .1 .3]]; 
S_col = repmat(S_col, ceil(K/size(S_col,1)),1); 
S_col = S_col(1:K, :); 
