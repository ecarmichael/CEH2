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

for ii = vid.NumFrames:-5:1
    frames(:,:,ii) = imgaussfilt(im2gray(read(vid, ii)), 1); 
end

% Calculate the frame differences
df = diff(double(frames),1,3); 




% sample movie
figure(101)
cla
for ii = 1:600
   imagesc(df(:,:, ii))
   pause(1/vid.FrameRate)
    drawnow
    
    
end








%%
out.pos_r = pos_r; 
out.fvec = fvec;
out.f_bin = f_val;
out.t_bin = (0:t_bin:ceil(pos_r.tvec(end))-2);