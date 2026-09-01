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

cfg_def.max_tbin_s = 0.001;  %bin size in seconds
cfg_def.max_t= .25;  %window size in seonds


% get the auto/cross correlation for each cell pair
fprintf('\nCell #    ')
for ii = length(S.t):-1:1
    for jj = length(S.t):-1:1

        [mono.cff{ii, jj} mono.t_vec] = ccf(cfg, S.t{ii}, S.t{jj});



    end
    fprintf('\b\b\b\b\b%2d/%2d', ii, length(S.t)); 
end

fprintf('\b\b\b\b\b%2d/%2d  - done\n', 0, length(S.t)); 

%%  plots?

if plot_flag

    figure(898)
    clf
    subplot(2,2,1)
    

    subplot(2,2,3)
    this_data = cell2mat(mono.cff(ii,:))';
    imagesc(mono.t_vec,1:length(S.t), this_data./max(this_data,[], 2))
    xlim([-.06 .06])

end