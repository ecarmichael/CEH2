function mono_metrics = MS_mono_metrics(S)
%% MS_mono_metrics:
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
% EC 2026-08-26   initial version 
%
%
%
%% initialize

bin_s = 0.001;  %bin size in seconds
win = .5;  %window size in seonds


% get the auto/cross correlation for each cell pair
fprintf('\nCell #    ')
for ii = length(S.t):-1:1
    for jj = length(S.t):-1:1

        [mono.cff{ii, jj} mono.t_vec] = ccf([], S.t{ii}, S.t{jj});



    end
    fprintf('\b\b\b\b\b%2d/%2d', ii, length(S.t)); 
end

fprintf('\b\b\b\b\b%2d/%2d  - done\n', 0, length(S.t)); 
