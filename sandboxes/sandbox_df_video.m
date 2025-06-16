% % [X, Y, Z] = extract_DF_F(,neuron.A,neuron.C,S_or,K_m+1)
% 
% 
% 
% 
% obj = VideoReader('msCam1.avi');
% vid = read(obj);
% frames = obj.NumberOfFrames;
% for x = 1 : frames
%     imwrite(vid(:,:,:,x),'msCam1.tif', 'WriteMode','append');
% end


%% df/f 

cam1 = squeeze(read(ms.vidObj{1})); 
%% flicker removal from : https://stackoverflow.com/questions/24365037/fix-fluorescent-flickering-from-video-using-matlab
cam1_s = nan(size(cam1)); 
for iF = 1:ms.vidObj{1}.NumFrames-1
    imDiff = cam1(:,:,iF+1) - cam1(:,:,iF);
    imDiff = medfilt2(imDiff); 
    imDiffRowSum = sum(imDiff,2); 
    dct_coeff = dct(imDiffRowSum); 
    flick_coeef = max(dct_coeff); 
    invdct = repmat(idct(flick_coeef),size(imDiff, 1), size(imDiff, 2));
    fix_diff = im2double(imDiff) - invdct; 
    cam1_s(:,:,iF) = im2double(cam1(:,:,iF+1)) + fix_diff; 
end

%% detrend? 
cam1_dt = NaN(size(cam1_s)); 
for ii = 1:size(cam1_s,1);
    
    for jj = 1:size(cam1_s, 2)
        
        cam1_dt(ii, jj,:) = detrend(cam1_s(ii, jj,:), 1, 0, 'SamplePoints', 33);  
        
    end
end


%%
min_frame = nanmean(cam1(:,200:700,:),3);
% all_min =min(cam1(150:400, 200:650,:), [], 'all'); 
d_cam1 = nan(size(cam1(:,200:700, :))); 
colormap('gray');
for iF = 1:ms.vidObj{1}.NumFrames- 50
%     min_idx = cam1(150:400, 200:650,iF) < min_frame;
% min_frame = nanmean(cam1_s(:,200:700,iF:iF+30),3);
   d_cam1(:,:,iF) = (nanmedian(cam1(:,200:700,iF:iF+5),3) ./ min_frame) -1; 
%     d_cam1_t(min_idx) = all_min; 
%      = d_cam1_t;
% if median(d_cam1(:,:,iF), 'all') < 0
%     d_cam1(:,:,iF) = zeros(size(d_cam1(:,:,iF))); 
% end
    
%     imwrite(d_cam1(:,:,iF),'d_cam1','gif','DelayTime',1/ms.vidObj{1}.FrameRate,'WriteMode','append');

end

%%
d_cam1_norm =1- d_cam1 ./ (nanmean(d_cam1, 'all') +(nanstd(d_cam1, [], 'all')*3)); 
% max_cur = (nanmean(d_cam1_norm, 'all') +(nanstd(d_cam1_norm, [], 'all')*1));
% min_cur = (nanmean(d_cam1_norm, 'all') -(nanstd(d_cam1_norm, [], 'all')*3));

clear G
for iF = 1:ms.vidObj{1}.NumFrames- 50
% d_cam1_norm(d_cam1_norm(:, :, iF) > max_cur) = max_cur; 
% d_cam1_norm(d_cam1_norm(:, :, iF) < min_cur) = min_cur; 
% d_cam1(:,:,iF) = normalize(d_cam1(:,:,iF));% ./(mean(d_cam1(:,:,iF), 'all') + std(d_cam1(:,:,iF), [], 'all')*2);

imagesc(d_cam1(:,:,iF)); 
% imagesc(cam1_s(:,:,iF)); 

    drawnow
    G(iF) = getframe(gcf) ;
    
end

%% make a gif?
Fs  =ms.vidObj{1}.FrameRate; 
for n = 1:2:200%ms.vidObj{1}.NumFrames-50
          [imind,cm] = rgb2ind(G(n).cdata,255);

      if n == 1
          imwrite(imind,cm,'temp_df.gif','gif','DelayTime',1/Fs, 'Loopcount',inf);
      else
          imwrite(imind,cm,'temp_df.gif','gif','DelayTime',1/Fs,'WriteMode','append');
      end
end

%%

cam2 = VideoReader('msCam2_df.avi'); 

cam2_v = read(cam2); 

% cam2_v = squeeze(rgb2gray(read(cam2))); 

cam2_g = NaN(cam2.Height, cam2.Width, cam2.NumFrames); 
for ii = 1:cam2.NumFrames

    cam2_g(:,:,ii) = im2double(rgb2gray(cam2_v(:,:,:,ii))); 
    
end
%%


min_frame = nanmean(cam2_g(:,:,:),3);
d_cam1 = nan(size(cam2_g(:,:, :))); 
colormap('gray');
for iF = 1:ms.vidObj{1}.NumFrames- 50
%     min_idx = cam1(150:400, 200:650,iF) < min_frame;
% min_frame = nanmean(cam1_s(:,200:700,iF:iF+30),3);
%    d_cam1(:,:,iF) = (nanmean(cam2_g(:,200:700,iF:iF+10),3) ./ min_frame) -1; 
%     d_cam1_t(min_idx) = all_min; 
d_cam1(:,:,iF) = medfilt2(cam2_g(:,:,iF)); 

%      = d_cam1_t;
% if median(d_cam1(:,:,iF), 'all') < 0
%     d_cam1(:,:,iF) = zeros(size(d_cam1(:,:,iF))); 
% end
    
%     imwrite(d_cam1(:,:,iF),'d_cam1','gif','DelayTime',1/ms.vidObj{1}.FrameRate,'WriteMode','append');

end
%%

    
%     d_cam1_norm =d_cam1 ./ (nanmean(d_cam1, 'all') +(nanstd(d_cam1, [], 'all')*3)); 
% max_cur = (nanmean(d_cam1_norm, 'all') +(nanstd(d_cam1_norm, [], 'all')*1));
% min_cur = (nanmean(d_cam1_norm, 'all') -(nanstd(d_cam1_norm, [], 'all')*3));
% d_min = 

for iF = 1:ms.vidObj{1}.NumFrames- 50
% d_cam1_norm(d_cam1_norm(:, :, iF) > max_cur) = max_cur; 
% d_cam1_norm(d_cam1_norm(:, :, iF) < min_cur) = min_cur; 


imagesc(d_cam1(:,:,iF)); 
% imagesc(cam1_s(:,:,iF)); 

    drawnow
   G(iF) =  getframe(gcf);
end

%% get data from full ms
idx = 13; 
rm_idx = 1:18;
rm_idx(idx) = []; 

    ms_temp = MS_remove_data_sandbox([], ms_seg_resize, rm_idx);

    ms_temp = MS_de_cell(ms_temp); 

%% try something else. 

cont_ms = ms_temp.SFPs(:,:,1);
c_ord = linspecer(size(ms_temp.SFPs,3)-200); 

figure(1)

norm_Raw = ms.RawTraces./norm(ms_temp.RawTraces,inf); 

for iC = size(ms_temp.RawTraces,2):-1:1
    norm_Raw(:,iC) = smooth( norm_Raw(:,iC), 15, 'sgolay');

%     norm_Raw(:,iC) = MS_norm_range(norm_Raw(:,iC), 0, 1);
end

norm_Raw(norm_Raw < 0) = 0; 

norm_Raw = norm_Raw ./max(norm_Raw, [], 'all'); 

% norm_Raw(norm_Raw >=1) = 1; 


norm_Raw(norm_Raw < 0.1) = 0; 



% norm_Raw = MS_norm_range(norm_Raw, 0, 1); 
clear F
F = NaN(1,length(round(size(ms_temp.RawTraces,1)/8))); 
parfor ii = 1:round(size(ms_temp.RawTraces,1)/8)
    disp(ii)
    
    rel_df = norm_Raw(ii,:);
    
    clf
    hold on
    
    for iC = size(ms_temp.SFPs,3)-200:-1:1
        [M] =contourc(ms_temp.SFPs(:,:,iC), [2,2]);
        [M2] =contourc(ms_temp.SFPs(:,:,iC), [4,4]);
%         clear c
        % c.ZData(c.ZData > 0) = iC;
        % cp = get(c, 'Children');
        % get(cp, 'cdata')
        % colormap(c_ord
        M(:,(M(1,:) < 5 | M(2,:) < 5)) = nan;
        M(1,:) = fillmissing(M(1,:), 'nearest');
        M(2,:) = fillmissing(M(2,:), 'nearest');
        
        M2(:,(M2(1,:) < 5 | M2(2,:) < 5)) = nan;
        M2(1,:) = fillmissing(M2(1,:), 'nearest');
        M2(2,:) = fillmissing(M2(2,:), 'nearest');
        if sum(rel_df) ==0
            continue
        end
        patch(M(1,:), M(2,:), c_ord(iC,:),'edgealpha', 0, 'FaceAlpha', rel_df(iC)/2)
        patch(M2(1,:), M2(2,:), c_ord(iC,:),'edgealpha', 0, 'FaceAlpha', rel_df(iC))

        % pause(.5)
    end
    
    xlim([0 ms_temp.width/ms_temp.ds]);
    ylim([0 ms_temp.height/ms_temp.ds]);
    set(gcf, 'color', 'k'); axis off; box on;

    
    drawnow;
    F(ii) = getframe(gcf) ;
    % plot(M(1,:), M(2,:), 'k')
end

%%
writerObj = VideoWriter(['temp_cell_color.avi']);
writerObj.FrameRate = mode(diff(ms.time));
writerObj.Quality = 100;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for iF = 1:length(F)
    % for ii =3962:2:5282
    
    % convert the image to a frame
    frame = F(iF) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

%% gif instead

Fs  =mode(diff(ms_temp.time)); 
for iF = 1:length(F)
      if iF == 1
          imwrite(F(iF).cdata,'temp_cell_df.gif','gif','DelayTime',1/Fs, 'Loopcount',inf);
      else
          imwrite(F(iF).cdata,'temp_cell_df.gif','gif','DelayTime',1/Fs,'WriteMode','append');
      end
end

%% manual video with df

v_id = 'msvideo.avi';

v_obj = VideoReader(v_id);

% v = read(v_obj)
%%
figure(19)
colormap('bone')

figure(20)
colormap('bone')
f_n = []; 

idx = floor(30*400):floor(30*420); 


mean_f = mean(squeeze(read(v_obj,[idx(1) idx(end)])),3);
max_f = mean(squeeze(read(v_obj,[idx(1) idx(end)])),'all');

d_f = NaN([size(mean_f), size(idx,1)]);  
d_f_s = d_f;  
F = cell(size(idx)); 
F2 = cell(size(idx)); 
F3 = cell(size(idx)); 

% mean_f = []; 
% for ii = floor(30*20):-1:floor(30*5)
%     
%     mean_f(:,:,ii) = mean(squeeze(read(v_obj,[ii-120 ii-1])),3);
% end

for ii = 1:length(idx)
    
%     mean_f = mean(squeeze(read(v_obj,[ii-120 ii-1])),3);

    this_f=  double(read(v_obj, idx(ii))) - mean_f; 
    

this_f = this_f ./ max_f; 
% this_f(this_f < .1) = 0; 
d_f(:,:,ii) = this_f;

% d_f(:,:,ii) = imgaussfilt(this_f, .2);

end

d_f_s = imgaussfilt3(d_f, .5); 

for ii = 1:length(idx)

% figure(19)
subplot(1,2,1)
r_f = double(read(v_obj, idx(ii))); 
    imagesc(r_f(crop(1):crop(2),crop(3):crop(4))); 
%     caxis([0 .5])
    
    axis off
%     F{ii} = getframe(gca) ;    

    subplot(1,2,2)
% figure(20)
    imagesc(d_f_s(crop(1):crop(2),crop(3):crop(4),ii)); 
    
%     montage({r_f(crop(1):crop(2),crop(3):crop(4)), d_f_s(crop(1):crop(2),crop(3):crop(4),ii)}); 
        caxis([0 .2])

   Square_subplots
axis off
    drawnow;
% 
% pause(.03)
    F3{ii} = getframe(gca) ;  
    
    
end

%%
Fs  =30; 
for iF = 1:length(F)
      if iF == 1
          imwrite(rgb2gray(F3{iF}.cdata),'temp_cell.gif','gif','DelayTime',1/Fs, 'Loopcount',inf);
      else
          imwrite(rgb2gray(F3{iF}.cdata),'temp_cell.gif','gif','DelayTime',1/Fs,'WriteMode','append');
      end
      
%       if iF == 1
%           imwrite(rgb2gray(F{iF}.cdata),'temp_cell_raw.gif','gif','DelayTime',1/Fs, 'Loopcount',inf);
%       else
%           imwrite(rgb2gray(F{iF}.cdata),'temp_cell_raw.gif','gif','DelayTime',1/Fs,'WriteMode','append');
%       end
end

%%

for ii = floor(30*20):-1:floor(30*5)
    
    f_n(:,:,ii) = double(read(v_obj, ii)); 
    
    f_df  = f_n(crop(1):crop(2),crop(3):crop(4),ii) - mean_f(crop(1):crop(2),crop(3):crop(4)); 
%     this_f(this_f < -2) = 0; 

% this_f = this_f ./ max_f; 
    
end
% mean_f = mean(f_n, 3); 

% f_df = f_n - mean_f; 

f_smooth = imgaussfilt3(f_n, .2); 

% f_df_smooth= imgaussfilt3(f_df, .2); 

for   ii = floor(30*20):-1:floor(30*5)
    
subplot(2,2,1)
% % v_img = double(read(v_obj, ii)); 
    imagesc(f_n(crop(1):crop(2),crop(3):crop(4),ii)); 
%     caxis([0 .5])
    
    subplot(2,2,2)
    imagesc(f_smooth(crop(1):crop(2),crop(3):crop(4),ii)); 
%     caxis([0 .5])
    
    subplot(2,2,3)
    imagesc(f_smooth(crop(1):crop(2),crop(3):crop(4),ii)); 
    
    subplot(2,2,4)
    imagesc(f_df_smooth(crop(1):crop(2),crop(3):crop(4),ii)); 
    drawnow;
%     pause(.03)
    F{end+1} = getframe(gcf) ;    
    
end

%%
writerObj = VideoWriter(['JKA_05_TFC.avi']);
writerObj.FrameRate = 30;
writerObj.Quality = 100;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for iF = 1:length(f_n)
    % for ii =3962:2:5282
    
    % convert the image to a frame
    frame = f_n(:,:,iF) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);



