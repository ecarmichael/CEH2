%% Example spatial plots using II data. 

cd('/home/williamslab/Dropbox (Williams Lab)/Williams Lab Team Folder/Eric/II_inter')

load('ck2cre_1377hd_2021_04_11_14_49_03_Open_Field.mat')


%% plot some simple things

iC = 1; % which cell to plot



figure(101)
subplot(2,2,1)
x_vals = data.SI(iC).cfg.X_bins; 
y_vals = data.SI(iC).cfg.Y_bins; 

imagesc(x_vals, y_vals, data.SI(iC).spatial.place.occupany)
title('Occupancy')
colorbar
subplot(2,2,2)
imagesc(x_vals, y_vals, data.SI(iC).spatial.place.posterior)
title('posterior prob.')
colorbar
subplot(2,2,3)
imagesc(x_vals, y_vals, data.SI(iC).spatial.place.Place_map_pval)
title('pval map')
colorbar
subplot(2,2,4)
imagesc(x_vals, y_vals, data.SI(iC).spatial.place.Sig_map)
title('Sig bins (scaled to max value in TC')
colorbar

