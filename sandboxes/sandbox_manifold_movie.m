%% manifold movie sandbox

% cd('/Volumes/Skadi/Using whole interpolate frames')

load('position_per_frame.mat');
load('v2.mat');
load('velocity_per_frame.mat');


%% plot an example
figure(1)
subplot(4,1,1:3)
scatter3(v2(:,2), v2(:,3), v2(:,4), 6*ones(1,length(position_per_frame)),position_per_frame, 'filled')

[az, el] = view;

this_dur = (33*60*2); % last value is number of minutes to run. 

F(this_dur) = getframe(gcf); % pre allocate; 

%%
% scatter3(v2(:,2), v2(:,3), v2(:,4), 6*ones(1,length(position_per_frame)),position_per_frame, 'filled')
c_ord = parula(ceil(max(position_per_frame))+1);
tvec = 1:length(position_per_frame); 
tvec = tvec*(1/30); 

figure('Position',[680 485 531 493],'MenuBar','none','ToolBar','none','resize','off')
set(gcf, 'position', [440   185   560   613])
for ii = 4:(30*60*2)+4
    subplot(4,1,1:3)
        hold on
%     scatter3(v2(:,2), v2(:,3), v2(:,4), 6*ones(1,length(position_per_frame)),position_per_frame)
%     alpha(.2)
    c1 = plot3(v2(ii-3:ii,2), v2(ii-3:ii,3), v2(ii-3:ii,4),'color', c_ord(round(position_per_frame(ii)+1),:), 'linewidth', 4);

    c3 = plot3(v2(ii,2), v2(ii,3), v2(ii,4), '.', 'color',[c_ord(ceil(position_per_frame(ii))+1,:) 0.25], 'markersize', 10);
    
    hold off
    xlim([min(v2(:,2)) max(v2(:,2))]);
    ylim([min(v2(:,3)) max(v2(:,3))]);
    zlim([min(v2(:,4)) max(v2(:,4))]);
    view(-39.394802263855850, 23.550427232944173);

    subplot(4,1,4)
    plot(position_per_frame(ii),0,'.', 'color', c_ord(round(position_per_frame(ii))+1,:), 'markersize', 30)
    hold on
    plot(position_per_frame(ii-2:ii),0, 'color', c_ord(round(position_per_frame(ii))+1,:), 'markersize', 30, 'linewidth', 4)
    hold off
    % add time. 
    t1 = text(10, .5, [num2str(tvec(ii),3) 's']);
        
    xlim([min(position_per_frame) max(position_per_frame)]);
    ylim([-1 1]); set(gca, 'ytick', []); 
    xlabel('position (cm)')
    set(gca, 'XDir','reverse')

     F(ii) = getframe(gcf) ;
      drawnow
      
      delete(c1);
      delete(t1); 
end

%% write to file. 
  writerObj = VideoWriter('Space_time_short.avi');
  writerObj.FrameRate = round(33);
  writerObj.Quality = 10; 
  % set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for ii=4:(30*60*2)+4
    % convert the image to a frame
    frame = F(ii) ;    
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);


%% manifold movie with Seq 