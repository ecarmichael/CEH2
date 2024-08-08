function spd = Spd_hist(velo, data_in,nShuff,  plot_flag)
%% compute the speed score using the Kropff method of pearson's correaltion relative to a shuffle. 

if nargin < 3
    nShuff = 100; 
    plot_flag = 0;
elseif nargin < 4
    plot_flag = 0;
end

nan_idx = isnan(velo) | isnan(data_in); 
velo(nan_idx) = []; 
data_in(nan_idx) = []; 


spd.R_corr = corr(velo, data_in, 'rows','pairwise', 'Type','Pearson');


% shuffle

shift_30 = 10;

spd.spd_shuff = NaN(1,nShuff); 

shuff_shifts = randperm(length(data_in) - shift_30, nShuff); 

for iS = nShuff:-1:1

    spd.spd_shuff(iS) = corr(data_in, circshift(velo, shuff_shifts(iS)+shift_30), 'Rows','pairwise', 'Type','Pearson');
end

spd.spd_99 = prctile(spd.spd_shuff, 99);
spd.spd_1 = prctile(spd.spd_shuff, 1);

spd.spd_mod = (spd.R_corr > spd.spd_99) || (spd.R_corr < spd.spd_1); 

% zscore
spd.spd_z = (spd.R_corr - mean(spd.spd_shuff))./(std(spd.spd_shuff));

% pval 
spd.spd_pval = sum(spd.spd_shuff > spd.R_corr,2)/nShuff;

%% Hass lab methods

% binned
spd.vd = (0:2.5:35);
spd.sp = [0, prctile(velo, 95)];

lincoeff  = polyfit(velo, data_in, 1); 

spd.y_m = polyval(lincoeff, velo, 2);

model_v_observed(:,1) = data_in; 
model_v_observed(:,2) = spd.y_m; 

[r, p] = corrcoef(model_v_observed);
spd.R = r(1,2); 
spd.P = p(1,2); 
spd.b = lincoeff(1); 
spd.y_int = lincoeff(2); 
spd.ym = spd.sp * spd.b + spd.y_int;

if plot_flag == 1
    figure(999)
    clf
    subplot(2,1,1)
    cla
    plot(velo, 'k')
    yyaxis right
    plot( data_in, 'r')
    
    subplot(2,1,2)  % add in smaller velo steps as well. 
    cla
    plot([0 35], [0 35]*(spd.b)+spd.y_int)
    hold on
    % plot in 1cm/s bins
    this_bin= 0:1:35; 
        [~,~,wb] = histcounts(velo, this_bin);
    for k = 1:max(this_bin)
        this_mv(k) = mean(data_in(wb==k));
    end
    plot(this_bin(1:end-1), this_mv,'ko','MarkerFaceColor','k', 'MarkerSize', 2.5)

    % plot in actual bins. 
    [~,~,wb] = histcounts(velo, spd.vd);
    for k = 1:max(wb)
        spd.mv(k) = mean(data_in(wb==k));
    end
    plot(spd.vd(1:end-1), spd.mv,'ko','MarkerFaceColor','r', 'MarkerSize', 10)

    xlabel('Running Speed (cm/s)');
    ylabel('Instantaneous Firing Rate (Hz)')
end

%% get the saturating R^2 value using function in CMBHOME.Utils.SaturatingRegression. 
f = fittype('d - a*exp(-b*x)','independent','x','coefficients',{'d' 'a' 'b'});
options = fitoptions(f);
options.Lower = [0 0 0];
options.Upper = [Inf,Inf,Inf];

[spd.sat.fit_params, goodness] = fit(data_in,velo,f,options);
spd.sat.rsquare = goodness.rsquare;

spd.sat.y_model = spd.sat.fit_params.d - spd.sat.fit_params.a*exp(-spd.sat.fit_params.b*data_in);


%% bonus plot
% 
% sub = 1; 
% epochs = table2array(unique(ArchtInhbSpdPwrTbl(ArchtInhbSpdPwrTbl.Subject==sub, 2))); 
% 
% All_R = NaN(size(epochs)); 
% All_spd_mod = All_R; 
% 
% for ii = length(epochs):-1:1
%     
%     this_t_amp = table2array(ArchtInhbSpdPwrTbl(ArchtInhbSpdPwrTbl.Subject==sub & ArchtInhbSpdPwrTbl.Epoch == epochs(ii), 3)); 
%     this_spd = table2array(ArchtInhbSpdPwrTbl(ArchtInhbSpdPwrTbl.Subject==sub & ArchtInhbSpdPwrTbl.Epoch == epochs(ii), 12)); 
% 
%     spd_out = Spd_hist(this_spd, this_t_amp, 100); 
%     
%     All_R(ii) = spd_out.R_corr; 
%     All_spd_mod(ii) = spd_out.spd_mod; 
%     
%     
%     
% end

% scatter(epochs(All_spd_mod==1), All_R(All_spd_mod==1), 'filled')
% hold on
% scatter(epochs(All_spd_mod==0), All_R(All_spd_mod==0))