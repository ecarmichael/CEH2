function norm_data = MS_norm_range(data, min_val, max_val)
%% MS_norm_range: normalize to range of values (min_val - max_val); 
%
%
%
%  EC initial version 2020
%
% based off of bz_NormToRange.m
% https://github.com/buzsakilab/buzcode/blob/05a8bd55ac48a6d50069c3f392d55fdc6e8cd5ec/visualization/bz_NormToRange.m

datarange = diff([min(data) max(data)]);
range= diff([min_val max_val]);

norm_data = (data-min(data))./datarange;
norm_data = norm_data.*range+min_val;
