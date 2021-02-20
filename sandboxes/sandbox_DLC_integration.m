% DLC integration sandbox



close all
restoredefaultpath
global PARAMS 

addpath(genpath('/home/ecarmichael/Documents/GitHub/CEH2')); 
%%
parent_dir = '/home/ecarmichael/Dropbox (Williams Lab)/Williams Lab Team Folder/Ingrid/Behav test and scripts/ck2cre-1359hd/2021_01_30/14_18_06'; 
% parent_dir = ('/mnt/Data/Behav test and scripts/ck-1361/2021_02_04'); 
cd(parent_dir); 
cd('BehavCam_1/')

%% load 

tracklength = 30;   % in cm  2021-01-30 is smallOF-30cm, 2021-02-02 is largeOF -60 sm

DLC_data = MS_collect_DLC(); 
DLC_fields = fieldnames(DLC_data);

tvec = table2array(readtable('timeStamps.csv', 'Headerlines', 1)); 

% get some video properties. 
vid_names = dir('*.avi');
temp= VideoReader(vid_names(1).name);
behav.height = temp.Height;
behav.width = temp.Width;

behav.dirName = parent_dir; 
behav.time = tvec(:,2); 
behav.fnum = tvec(:,1); 
behav.buffer = tvec(:,3); 
behav.Fs = mode(diff(behav.time)); 

% main position using LED
position = DLC_data.LED(:,1:2); 

% scale pixels to cm. 
    position = position*tracklength/behav.width;

    % this section is straight out of msExtractBehavior by GE; 
    time = behav.time(~isnan(position(:,1)));
    position = position(~isnan(position(:,1)),:);
    [time, index] = unique(time); % GE added
    position_q(:,1) = interp1(time,position(index,1),behav.time); %GE renamed from: position = interp1(time,position,behav.time);
    position_q(:,2) = interp1(time,position(index,2),behav.time); % Added by GE
    dt = median(diff(behav.time/1000)); 
    position_q(:,1) = smooth(position_q(:,1),ceil(1/dt)); % GE renamed from: position = smoothts(position','b',ceil(1/dt))';
    position_q(:,2) = smooth(position_q(:,2),ceil(1/dt)); % Added by GE
    behav.position = position_q;
    
    dx = [0; diff(position_q(:,1))]; % GE modified
    dy = [0; diff(position_q(:,2))]; % GE modified

    behav.speed = sqrt((dx).^2+(dy).^2)/dt;
    behav.speed = smooth(behav.speed,ceil(1/dt));
    behav.speed(1) = behav.speed(2); % fill in the first with the next point as an estimate. 
    
    % same smoothing for other DLC positions
    
    for iF = 1:length(DLC_fields)
        behav.(DLC_fields{iF}) = [];
        this_pos = [];
        
        this_pos = DLC_data.(DLC_fields{iF})(:,1:2)*tracklength/behav.width; % pixel -> cm correction.
        
        this_pos_q(:,1) = interp1(time,this_pos(index,1),behav.time); %GE renamed from: position = interp1(time,position,behav.time);
        this_pos_q(:,2) = interp1(time,this_pos(index,2),behav.time); % Added by GE
        this_pos_q(:,1) = smooth(this_pos_q(:,1),ceil(1/dt)); % GE renamed from: position = smoothts(position','b',ceil(1/dt))';
        this_pos_q(:,2) = smooth(this_pos_q(:,2),ceil(1/dt)); % Added by GE
        
        behav.(DLC_fields{iF}) =  this_pos_q; 
        
    end

% get the HD if there are ear trackers
if sum(contains(lower(DLC_fields), 'nose'))==1 && sum(contains(lower(DLC_fields), 'r_ear'))==1  && sum(contains(lower(DLC_fields), 'l_ear')) == 1
    fprintf('<strong>Found ear tracking in DLC files.  Extracting HD...</strong>\n')
   
        % compute midpoint between ears
        ear_mid(:,1) = (behav.R_Ear(:,1) + behav.L_Ear(:,1))/2; % get x mid
        ear_mid(:,2) = (behav.R_Ear(:,2) + behav.L_Ear(:,2))/2; % get x mid
       
        behav.HD = rad2deg(atan2(ear_mid(:,2) - behav.LED(:,2),ear_mid(:,1) - behav.LED(:,1))); 
        
        % angle test plot
%                
%         hold on
%          for iF =1:length(behav.Nose)
%         plot(behav.position(iF,1), behav.position(iF,2), 'or')
%         plot([behav.R_Ear(iF,1),  behav.L_Ear(iF,1)], [behav.R_Ear(iF,2),  behav.L_Ear(iF,2)], 'sr')
%         plot(ear_mid(iF, 1), ear_mid(iF,2), '.b')
%         plot(behav.Nose(iF,1), behav.Nose(iF,2), 'xg')
%         plot([ear_mid(iF,1), behav.LED(iF,1)],[ear_mid(iF,2), behav.LED(iF,2)], 'k')
%         xlim([min(behav.position(:,1)) max(behav.position(:,1))]); 
%         ylim([min(behav.position(:,2)) max(behav.position(:,2))]);
%         text(0, 0, num2str(behav.HD(iF)));
%         text(0, 5, num2str(behav.speed(iF))); 
%         drawnow
%         pause(0.0303/30)
%                 cla(gca)
%          end
end



        
% test plot
% figure(111)

% plot(behav.position(:,1), behav.position(:,2), '.')

% write back the file
save([parent_dir filesep 'behav_DLC.mat'], 'behav')
fprintf('Output saved as ''behav_DLC.mat'' in %s\n', parent_dir)


