function hd_temp = MS_Load_v4HD(data_dir, plot_flag)
%% MS_Load_v4HD:
%
%
%
%    Inputs: 
%    - data_dir: [string] path to data
%
%
%
%    Outputs: 
%    - hd: [struct] pith roll and yaw in the tsd format. 
%
%
%
%
% EC 2022-12-19   initial version 
%
%
%
%% initialize

if nargin < 1
    data_dir = cd;
end

cd(data_dir);

%% load the data
hd_temp = [];

hd_tbl = readtable('headOrientation.csv');
hd_temp.time = hd_tbl.TimeStamp_ms_ ./ 1000; 
hd_temp.quat = [hd_tbl.qw, hd_tbl.qx, hd_tbl.qy, hd_tbl.qz]; 

% convert to euclidean 
eulZYX = quat2eul(hd_temp.quat);


eulZYX_df_t(:,1) = conv2(diff(eulZYX(:,1)),gausswin(60, 3),'same'); 
eulZYX_df_t(:,2) = conv2(diff(eulZYX(:,2)),gausswin(60, 3),'same');
eulZYX_df_t(:,3) = conv2(diff(eulZYX(:,3)),gausswin(60, 3),'same');

eulZYX_df = [eulZYX_df_t; eulZYX_df_t(end,:)];

hd_temp.motion = sqrt(movmean(sum(eulZYX_df,2).^ 2, 5)); 
%% plot

if plot_flag
    
   figure(919)
   subplot(3,1,1)
   plot(hd_temp.time, eulZYX)
   ylabel('Euclidean movement')
   legend({'Yaw (Z)', 'Pitch (Y)', 'Roll (X)'}); 
   
   subplot(3,1,2)
      plot(hd_temp.time, hd_temp.motion)
ylabel('Movement (sqrt movmean: 5 samples)')

   subplot(3,1,3)
   plot(hd_temp.time, hd_temp.quat)
       ylabel('Quaternion movement')
       xlabel('time (s)')

    
end

%% export in tsd format

hd = tsd(hd_temp.time, [eulZYX, hd_temp.motion]', {'Z', 'Y', 'X', 'Motion'}); 
