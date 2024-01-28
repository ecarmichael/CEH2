function[z_x_c, z_ac_max, mean_x_c, ac_max] = MS_Asmbly_xcor(Proj_in, tvec, R_thresh, xc_bin, t_max)

if nargin < 3
    R_thresh = 5;
    xc_bin = 0.1; in seconds
    t_max = 2.5;
elseif nargin <4
    xc_bin = 0.1; in seconds
    t_max = 2.5;
elseif nargin < 5
    t_max = 2.5;
end

cfg.binsize = xc_bin;
cfg.max_t = t_max;

%% convert the data to timestamps; 

for ii = size(Proj_in, 1):-1:1
    [~, this_ts] = findpeaks(Proj_in(ii,:), 'MinPeakHeight', R_thresh); 
    
    TS{ii} = tvec(this_ts); 
end

%% create a null distibution


for iS = 1000:-1:1
    idx = randsample(1:size(Proj_in,1),2, false ); 
    this_proj = Proj_in(idx,:);

    this_shuff = circshift(this_proj(2,:), randsample(1:length(Proj_in), 1));
    
    [s_ac, lag] = xcov(this_shuff, this_proj(1,:),cfg.max_t/mode(diff(tvec)),  'coeff');
    s_ac_max(iS) = max(s_ac);
    
    TS_shuff = randsample(tvec, length(TS{idx(1)}));
    TS_shuff2 = randsample(tvec, length(TS{idx(2)}));

    shuff_x_c(iS) = mean(ccf(cfg, TS_shuff, TS_shuff2));
    
end
        

%% get the xcorr


for ii = size(Proj_in, 1):-1:1
    
    for jj = size(Proj_in, 1):-1:1
        
        % create a shuffle null 

        
        [x_c{ii, jj}, t_vec] = ccf(cfg, TS{ii}, TS{jj});
        mean_x_c(ii, jj) = mean(x_c{ii,jj}); % as per wilson mcnaughton; 
        z_x_c(ii, jj) = (mean_x_c(ii,jj) - mean(shuff_x_c)) / std(shuff_x_c);
        
        [ac, lag] = xcov(Proj_in(ii,:), Proj_in(jj,:),cfg.max_t/mode(diff(tvec)),  'coeff');
        if nansum(ac) == 0
            ac_max(ii, jj) = NaN;
            lag_max(ii, jj) = NaN;
        else
            [ac_max(ii,jj),idx] = max(ac);
            lag_max(ii,jj) = lag(idx);
        end
        z_ac_max(ii,jj) = (mean(ac_max(ii,jj)) - mean(s_ac_max)) / std(s_ac_max); 
        
    end
end


%%
