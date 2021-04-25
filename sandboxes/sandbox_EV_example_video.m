addpath(genpath('/home/ecarmichael/Documents/GitHub/CEH2')); 


inter_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\10.Manifold';
addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));
addpath('C:\Users\ecarm\Documents\GitHub\OASIS_matlab')
oasis_setup;
disp('OASIS added for deconv');

raw_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\RAW Calcium\LT&sleep\12_15_2019_540day5\H11_M52_S58lt';
ms_dir = 'C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eva\final_analysis\540\12_15_2019_540day5'; 

%% 
clearvars -except raw_dir ms_dir
cd(raw_dir)


b_obj = VideoReader('behavCam2.avi');
for iF = b_obj.NumFrames:-1:1
    b_f{iF} = read(b_obj, iF);
    b_t(iF) = b_obj.CurrentTime; 
end

ca_obj = VideoReader('msCam2.avi');
for iF = ca_obj.NumFrames:-1:1
    ca_f{iF} = read(ca_obj, iF);
    ca_t(iF) = ca_obj.CurrentTime; 
end

%%
cd(ms_dir)
load('ms.mat')
ms = msExtractBinary_detrendTraces(ms); 
frame_offset = ms.timestamps(1); 
Fs = mode(diff(ms.time)); 
%% deconv debugger. 
% % counter_init(size(ms.RawTraces,2),size(ms.RawTraces,2))
% for iChan = 25:-1:1
%     counter(iChan, size(ms.RawTraces,2))
%     [denoise,deconv] = deconvolveCa(ms.detrendRaw(:,iChan), 'foopsi', 'ar1', 'smin', -3, 'optimize_pars', true, 'optimize_b', true);
%     ms.denoise(:,iChan) = denoise; 
%     ms.deconv(:,iChan) = deconv; 
% end

%% debugging
if ishandle(111)
    close(111)
end
figure(111)
% iChan = 25;
% hold on
% plot(ms.detrendRaw(:,iChan), 'k');
% plot(ms.denoise(:,iChan)-.2, 'r');
% plot(((ms.deconv(:,iChan)./max(ms.deconv(:,iChan)))*.1) -.2, 'b');
% plot((ms.Binary(:,iChan)*.1)+max(ms.RawTraces(:,iChan)), 'g'); 
% MS_plot_ca_trace(ms.FiltTraces(1:33*60,1:50)')

%% plot stuff
close all
figure(101)
subplot(3,4,1:4)
imagesc(b_f{1})
title(sprintf('Behav Time = %.3fs | Ca Time = %.3fs', b_t(1), ca_t(1)));
set(gca, 'xtick', [], 'ytick', [])

subplot(3,4,[5 6 9 10])
MS_Ca_Raster(ms.Binary(1:ms.timestamps(2),:)'); %, ms.time(1:ms.timestamps(2))'/1000
set(gca, 'xtick', [])
ylabel('cell number')


subplot(3,4,[ 7 8 11 12])
imagesc(ca_f{1})
colormap('gray')
set(gca, 'xtick', [], 'ytick', [])

set(gcf, 'position', [511  43 1242 935],'MenuBar','none','ToolBar','none')

clear F
for iF = Fs:4:size(b_t,2)-Fs
    
    subplot(3,4,1:4)
    imagesc(b_f{iF})
    title(sprintf('Behav Time = %.3fs | Ca Time = %.3fs', b_t(iF), ca_t(iF)));

    subplot(3,4,[5 6 9 10])
    xlim([ms.timestamps(1)-Fs + iF  ms.timestamps(1)+Fs + iF])
    h = vline(median([ms.timestamps(1)-Fs + iF  ms.timestamps(1)+Fs + iF]));

    subplot(3,4,[ 7 8 11 12])
    imagesc(ca_f{iF})

    drawnow;
    F(iF) = getframe(gcf) ;
    delete(h);

end

%% limited cell version

cells = 1:8;
sub_mat = reshape(1:((2+length(cells))*4),4,2+length(cells))'; % matrix to pull subplot values from

figure(108)
subplot(2+length(cells),4,1:8)
imagesc(b_f{1})
title(sprintf('Behav Time = %.3fs | Ca Time = %.3fs', b_t(1), ca_t(1)));
set(gca, 'xtick', [], 'ytick', [])

subplot(3,4,[5 6 9 10])
% MS_Ca_Raster(ms.Binary(1:ms.timestamps(2),:)'); %, ms.time(1:ms.timestamps(2))'/1000
set(gca, 'xtick', [])
ylabel('cell number')


subplot(2+length(cells),4,[ 7 8 11 12])
imagesc(ca_f{1})
colormap('gray')
set(gca, 'xtick', [], 'ytick', [])

set(gcf, 'position', [511  43 1242 935],'MenuBar','none','ToolBar','none')


%%
writerObj = VideoWriter(['C:\Users\ecarm\Dropbox (Williams Lab)\Williams Lab Team Folder\Eric\EV_inter' filesep 'Example_Ca_behav.avi']);
writerObj.FrameRate = 16.5;
writerObj.Quality = 40;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for iF = Fs:4:size(b_t,2)-Fs
    % for ii =3962:2:5282
    
    % convert the image to a frame
    frame = F(iF) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


