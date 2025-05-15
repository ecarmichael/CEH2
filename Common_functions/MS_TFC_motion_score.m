function out = MS_TFC_motion_score(data_dir, threshold)
%% MS_TFC_motion_score:
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
% EC 2025-04-24   initial version 
%
%
%
%% initialize

% eztrack suggested a movement threshold of 2 x 99.99th percentile of the
% animal free arena. since we do't have this, ...

mo_thresh = 100; 

vid = VideoReader('behav.avi'); 

frames = zeros(vid.Height, vid.Width, vid.NumFrames, 'uint8');

for ii = vid.NumFrames:-1:1
    frames(:,:,ii) = imgaussfilt(im2gray(read(vid, ii)), 1); 
end

% Calculate the frame differences
df = diff(double(frames),1,3); 

m_df = []; keep_idx = false(size(df,3),1); 
for ii = 600:-2:1
    z = (mean(cat(3,abs(df(:,:,ii-1)), abs(df(:,:,ii))),3)); 
    m_df(ii) = mean(z, 'all'); 
    keep_idx(ii) = true;
end


%%
% sample movie
figure(101)
clf
len = 600; 
   subplot(4,1,4)
   hold on
plot(1:len, m_df(1:len))
for ii = 1:len
    subplot(4,1,1:3)

   imagesc(df(:,:, ii))
   pause(1/vid.FrameRate)
   
   subplot(4,1,4)
   delete(xl); 
   xl = xline(ii);
    drawnow
    
    
end

%% get a motion vector

% figure(






%%
out.pos_r = pos_r; 
out.fvec = fvec;
out.f_bin = f_val;
out.t_bin = (0:t_bin:ceil(pos_r.tvec(end))-2);