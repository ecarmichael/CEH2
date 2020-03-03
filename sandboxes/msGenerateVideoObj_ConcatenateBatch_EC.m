function ms = msGenerateVideoObj_ConcatenateBatch_EC(dirName, filePrefix)

    MAXFRAMESPERFILE = 1000; %This is set in the miniscope control software
    ms.dirName = dirName;
    % find avi and dat files    
    if ~isempty(strfind(filePrefix,'msCam'))
        aviFiles = dir(fullfile(dirName, '**', 'msCam*.avi')); %looks at whats inside this directory to see files
    elseif ~isempty(strfind(filePrefix,'behavCam'))
        aviFiles = dir(fullfile(dirName, '**', 'behavCam*.avi')); %looks at whats inside this directory to see files
    end 
    
    datFiles = dir(fullfile(dirName, '**', 'timestamp*.dat'));    
    
    %we need to sort the aviFiles and the datFiles such that they are in
    %the order that matches the current folder: 
    timeLocation = strfind(datFiles(1).folder,[filesep 'H']); % fixed for use across unix/pc
    s = strings(length(datFiles),1);
    folderNames = strings(length(datFiles),1);
    for folderNum = 1 : length(datFiles)        
        folderNames(folderNum,1)  = datFiles(folderNum).folder ; 
        m =folderNames{folderNum}(timeLocation+1:end); 
        s(folderNum,1) = m; 
    end 
    
    for num =1: length(datFiles)
        order{1,num} = sscanf(s(num,1),'H%d_M%d_S%d*');
    end 
    
      %this next chunck takes the order variable which has the hours,
      %minutes, and seconds saved, and create a new array with the format
      %of hours, minutes, seconds saved in a single array, in a format the
      %computer can properly read and sort. 
      orderTime = strings(length(datFiles),1);
      for num = 1: length(orderTime)
          %001
          if nnz(isstrprop(num2str(order{1,num}(1)),'digit')) > 1 && nnz(isstrprop(num2str(order{1,num}(2)),'digit')) > 1  && nnz(isstrprop(num2str(order{1,num}(3)),'digit'))== 1 
              orderTime(num,1) = sprintf('%d:%d:0%d',order{1,num}(1),order{1,num}(2),order{1,num}(3));
         %010    
          elseif nnz(isstrprop(num2str(order{1,num}(1)),'digit')) > 1 && nnz(isstrprop(num2str(order{1,num}(2)),'digit')) == 1  && nnz(isstrprop(num2str(order{1,num}(3)),'digit')) > 1 
              orderTime(num,1) = sprintf('%d:0%d:%d',order{1,num}(1),order{1,num}(2),order{1,num}(3)); 
          %011    
          elseif nnz(isstrprop(num2str(order{1,num}(1)),'digit')) > 1 && nnz(isstrprop(num2str(order{1,num}(2)),'digit')) == 1  && nnz(isstrprop(num2str(order{1,num}(3)),'digit')) == 1 
              orderTime(num,1) = sprintf('%d:0%d:0%d',order{1,num}(1),order{1,num}(2),order{1,num}(3));              
          %100    
          elseif nnz(isstrprop(num2str(order{1,num}(1)),'digit')) == 1 && nnz(isstrprop(num2str(order{1,num}(2)),'digit')) > 1  && nnz(isstrprop(num2str(order{1,num}(3)),'digit')) > 1 
              orderTime(num,1) = sprintf('0%d:%d:%d',order{1,num}(1),order{1,num}(2),order{1,num}(3));
           %101    
          elseif nnz(isstrprop(num2str(order{1,num}(1)),'digit')) == 1 && nnz(isstrprop(num2str(order{1,num}(2)),'digit')) > 1  && nnz(isstrprop(num2str(order{1,num}(3)),'digit'))== 1 
              orderTime(num,1) = sprintf('0%d:%d:0%d',order{1,num}(1),order{1,num}(2),order{1,num}(3)); 
          %110    
          elseif nnz(isstrprop(num2str(order{1,num}(1)),'digit')) == 1 && nnz(isstrprop(num2str(order{1,num}(2)),'digit')) == 1  && nnz(isstrprop(num2str(order{1,num}(3)),'digit')) > 1  
              orderTime(num,1) = sprintf('0%d:0%d:%d',order{1,num}(1),order{1,num}(2),order{1,num}(3));              
          %111   
          elseif nnz(isstrprop(num2str(order{1,num}(1)),'digit')) == 1 && nnz(isstrprop(num2str(order{1,num}(2)),'digit')) == 1  && nnz(isstrprop(num2str(order{1,num}(3)),'digit')) == 1 
              orderTime(num,1) = sprintf('0%d:0%d:0%d',order{1,num}(1),order{1,num}(2),order{1,num}(3));
          %000    
          else 
              orderTime(num,1) = sprintf('%d:%d:%d',order{1,num}(1),order{1,num}(2),order{1,num}(3)); 
      end 
      end      
     
     % now rearranging the original datFiles and aviFiles: 
     [~,idxFolder] = sort(orderTime);
     for num =1: length(datFiles)
         new_datFiles(num) =  datFiles(idxFolder(num)); 
         new_folderNames(num) = folderNames(idxFolder(num)); 
     end 
     datFiles = new_datFiles; 
     folderNames = new_folderNames'; 
     
     aviFolders = strings(length(aviFiles),1);
     for num =1 : length(aviFiles)
         aviFolders(num) = aviFiles(num).folder; %we can index easier this way 
     end 
     
     
     for num =1: length(folderNames)
        idx_aviFiles = find(aviFolders == folderNames(num)); 
        if num == 1 
            new_aviFiles = aviFiles(idx_aviFiles); 
        else
            new_aviFiles = [new_aviFiles; aviFiles(idx_aviFiles)];
        end 
     end 
       aviFiles = new_aviFiles; 
    
       
       %if the one before it is more than 1 away and not a 1, then move. 
       for videoNum =1: length(aviFiles)
        if ~isempty(strfind(filePrefix,'msCam'))
            videoOrder(videoNum) = sscanf(aviFiles(videoNum).name,'msCam%d.avi');
        elseif ~isempty(strfind(filePrefix,'behavCam'))
            videoOrder(videoNum) = sscanf(aviFiles(videoNum).name,'behavCam%d.avi');
        end           
       end
       videoOrder = videoOrder'; 
       
        for num =1 : length(aviFiles)
         aviFolders(num) = aviFiles(num).folder; %we can index easier this way 
        end 
     
       for num = 1: length(folderNames)
              location = strfind(aviFolders, folderNames(num)); %find the folder we are looking at 
              idx_location = find(~cellfun(@isempty,location)); %the index of the videos in that folder              
              videos = aviFiles(idx_location);  %save them in a separate array to sort and put back into aviFiles
              videos_idx = videoOrder(idx_location);
              
              [~,idxmatr] = sort(videos_idx); %sort the videos with their saved index
              videos = videos(idxmatr); %with their saved index reorder "videos"
              
              %Now placing them back into aviFiles:              
              aviFiles = [aviFiles(1:idx_location(1)-1);    videos;    aviFiles(idx_location(end)+1:end)];              
       end    
       
% need to make folder save the locations of the avi files
    %now we need to rename these so that these files are access correctly: 
    for ii =1: (length(aviFiles) +length(datFiles))  
        if ii > length(aviFiles)
            folder(ii) = datFiles(ii-length(aviFiles)); 
        else
            folder(ii) = aviFiles(ii); 
        end 
    end 
%we need to change the names of each of the video files:
current_path = cd; 
folderNames = dir; 
count1 = 1; 
count2 =length(aviFiles)+1;
    %gets the names of all the folder names in the directory       
            folderNames(1) = []; 
            folderNames(1) = [];
            tempNames = folder;        
     
    ms.numFiles = 0;        %Number of relevant .avi files in the folder
    ms.numFrames = 0;       %Number of frames within said videos
    ms.vidNum = [];         %Video index
    ms.frameNum = [];       %Frame number index
    ms.maxFramesPerFile = MAXFRAMESPERFILE; %finds the maximum number of frames contained in a single throughout all videos
    ms.dffframe = [];
    
    %find the total number of relevant video files
    ms.numFiles = length(aviFiles);   
    
    Timestampmin = 1000000;     %Arbitrary VERY large number for min value storage, must be larger than datenum stamp in folder
    Timestamp = 0;              %Total number of frames observed (including previous video itterations)
    o = NaN(1,ms.numFiles);     %index movie file order
    anomilynames = NaN(1,4);
    count = 0;
    
    o = nan(1,length(folder));
    stamps = zeros(1,length(folder));
    j = 1;
    k =1;
    for ii = 1:length(folder)
        if ~folder(ii).isdir && ~isempty(strfind(folder(ii).name,filePrefix))
            o(ii) = k;
            k = k+1;
        elseif ~folder(ii).isdir && ~isempty(strfind(folder(ii).name,'timestamp'))
            stamps(ii) = j;
            j=j+1;
        end
    end
    [a b] = sort(o);
    folder = folder(b);
    o = o(b);
    [a b] = sort(stamps);
    folder = folder(b);
    stamps = stamps(b);
    
    %generate a vidObj for each video file. Also calculate total frames
    for ii=1:ms.numFiles
        j = o(1,ii);                                                     %call on .avi files chronologically   
        ms.vidObj{ii} = VideoReader([folder(j).folder filesep folder(j).name]); %Read .avi video file
        ms.vidNum = [ms.vidNum ii*ones(1,ms.vidObj{ii}.NumberOfFrames)];  %Store video index into ms for future use outside this fn
        ms.frameNum = [ms.frameNum 1:ms.vidObj{ii}.NumberOfFrames];      %Current frame # in total
        ms.numFrames = ms.numFrames + ms.vidObj{ii}.NumberOfFrames;      %Total number of frames
    end
    ms.height = ms.vidObj{1}.Height;        %video dimentions
    ms.width = ms.vidObj{1}.Width;
    
    camNum = [];
    frameNum = [];
    sysClock = [];
    buffer1 = [];
    frameTot1 = 0;
    frameTot0 = 0;
    %% read timestamp information
    for ii=1:length(datFiles)
        idx1 = []; 
        idx0 = [];
        camNumtemp = []; 
        dataArray =[]; 
        fileID =[];
        
        if strfind(datFiles(ii).name,'timestamp')     
            if isunix
                fileID = fopen(sprintf('%s//%s',datFiles(ii).folder,datFiles(ii).name)); %looks for timestamp, but will not iterate through because it'll only accept the first timestamp
            else
                fileID = fopen(sprintf('%s\\%s',datFiles(ii).folder,datFiles(ii).name)); %looks for timestamp, but will not iterate through because it'll only accept the first timestamp
            end
            dataArray = textscan(fileID, '%f%f%f%f%[^\n\r]', 'Delimiter', '\t', 'EmptyValue' ,NaN,'HeaderLines' ,1, 'ReturnOnError', false);    %read file and make sure it is not empty
            camNumtemp = dataArray{:, 1};       %camera number
            frameNumtemp = dataArray{:, 2};            
            A = unique(camNumtemp); 
            sysClocktemp = dataArray{:, 3};     %system clock 
            buffer1temp = dataArray{:, 4};      %buffer
            
            % hold the timestamp variables in the ms output (saves time
            % later when trying to align MS to NLX
            array_temp = dataArray;
            array_temp(5) = [];
            ms.timestamp_data{ii} = cell2table(array_temp,'Variablenames',{ 'camNum', 'frameNum', 'sysClock', 'buffer1'});
            clear array_temp
            %--------------------------------------------------------------
            if ii ==1
                frameTot = 0; 
                frameTot1 = 0;
                frameTot0 = 0; 
                sysClock = [];
                sysClock1 =[];
                sysClock0 =[];
                sysClocktemp1 =[];
                sysClocktemp0 =[];                
            end 
            
            %--------------------------------------------------------------
            if length(A) == 1
                idx = find(camNumtemp == A(1));                      
                camNum1 = frameNumtemp(idx); 
                cam1 = max(idx); 
                clearvars dataArray;            %clear variables from dataArray
                fclose(fileID);
                sysClocktemp(1) = 0;  
                ms.timestamps0(ii,1) = max(frameNumtemp(find(camNumtemp == A(1)))); 
                ms.timestamps1(ii,1) = max(frameNumtemp(find(camNumtemp == A(1))));
                    if ii == 1
                        camNum = camNumtemp ;
                        %frameNum =  frameNumtemp;
                        sysClock1 = sysClocktemp(idx);
                        sysClock0 = sysClocktemp(idx);
                        sysClocktemp1 = sysClocktemp(cam1);
                        sysClocktemp0 = sysClocktemp(cam1);
                        buffer1 = buffer1temp; 
                        frameTot1 = frameTot1; 
                        frameTot0 = frameTot0; 
                    else
                        camNum = [camNum; camNumtemp] ;
                        frameNum = [frameNum; frameNumtemp];
                        sysClock1 = [sysClock1; (sysClocktemp(idx) + sysClocktemp1)]; %we dont yet add times together to get the end ms.time result. 
                        sysClock0 = [sysClock0; (sysClocktemp(idx) + sysClocktemp0)];                        
                        sysClocktemp1 = sysClocktemp(cam1) + sysClocktemp1;
                        sysClocktemp0 = sysClocktemp(cam1) + sysClocktemp0;
                        buffer1 = [buffer1; buffer1temp];
                    end
                        frameTot1 = frameTot1 + camNum1(length(camNum1));
                        frameTot0 = frameTot0 + camNum1(length(camNum1));
                
            %--------------------------------------------------------------
            else 
                sysClocktemp(1) = 0; 
                sysClocktemp(2) = 0; 
                idx0 = find(camNumtemp == min(A));     
                idx1 = find(camNumtemp == max(A));     
                camNum0 = frameNumtemp(idx0);
                camNum1 = frameNumtemp(idx1);                               
                cam1 = max(idx1);                              
                cam0 = max(idx0);   
%                 sysClocktempidx1 = sysClocktemp(idx1);
%                 sysClocktempidx0 = sysClocktemp(idx0);
                clearvars dataArray;            %clear variables from dataArray
                fclose(fileID);
                ms.timestamps0(ii,1) = max(frameNumtemp(find(camNumtemp == A(1)))); 
                ms.timestamps1(ii,1) = max(frameNumtemp(find(camNumtemp == A(2))));
                
                    if ii == 1
                        camNum = camNumtemp ;                         
                        %frameNum =  frameNumtemp;
                        sysClock1 = sysClocktemp(idx1);
                        sysClock0 = sysClocktemp(idx0);
                        buffer1 = buffer1temp;
                        sysClocktemp1 = sysClocktemp(cam1);
                        sysClocktemp0 = sysClocktemp(cam0);
                        frameTot1 = frameTot;
                        frameTot0 = frameTot;
                    else
                        camNum = [camNum; camNumtemp] ;
                        frameNum = [frameNum; frameNumtemp];                       
                        buffer1 = [buffer1; buffer1temp];
                        sysClock1 = [sysClock1; (sysClocktemp(idx1)+ sysClocktemp1)];
                        sysClock0 = [sysClock0; (sysClocktemp(idx0)+ sysClocktemp0)];
                        sysClocktemp1 = sysClocktemp(cam1) + sysClocktemp1;
                        sysClocktemp0 = sysClocktemp(cam0) + sysClocktemp0;
                    end
                    frameTot1 = frameTot1 + camNum1(length(camNum1));
                    frameTot0 = frameTot0 + camNum0(length(camNum0)); 
            end              
        end     
    end 
        
    
    for j=max(camNum):-1:0   
        if (sum(camNum==j)~=0)
            if (frameTot1 == ms.numFrames) && (sum(camNum==j) == length(camNum(camNum==j)));
                ms.camNumber = j;
                ms.time = sysClock1;
                ms.time(1) = 0;
                ms.maxBufferUsed = max(buffer1(camNum==j));
                ms.timestamps = ms.timestamps1; 
            elseif (frameTot0 == ms.numFrames) && (sum(camNum==j) == length(camNum(camNum==j)));
                ms.camNumber = j;
                ms.time = sysClock0;
                ms.time(1) = 0;
                ms.maxBufferUsed = max(buffer1(camNum==j));
                ms.timestamps = ms.timestamps0; 
            else
                display(['Problem matching up timestamps for ' dirName]);
            end
        end
        
        
         if strcmp(datFiles(ii).name, 'settings_and_notes.dat')
            fileID = fopen([dirName filesep datFiles(ii).name],'r');
            textscan(fileID, '%[^\n\r]', 1, 'ReturnOnError', false);
            dataArray = textscan(fileID, '%s%s%s%s%[^\n\r]', 1, 'Delimiter', '\t', 'ReturnOnError', false);
            ms.Experiment = dataArray(:,1);
            ms.Experiment = string(ms.Experiment{1});
        end
    end

    
%     %figure out date and time of recording if that information if available
%     %in folder path
    idx = strfind(dirName, '_');
    idx2 = strfind(dirName,filesep);
    if (length(idx) >= 4)
        ms.dateNum = datenum(str2double(dirName((idx(end-2)+1):(idx2(end)-1))), ... %year
            str2double(dirName((idx2(end-1)+1):(idx(end-3)-1))), ... %month
            str2double(dirName((idx(end-3)+1):(idx(end-2)-1))), ... %day
            str2double(dirName((idx2(end)+2):(idx(end-1)-1))), ...%hour
            str2double(dirName((idx(end-1)+2):(idx(end)-1))), ...%minute
            str2double(dirName((idx(end)+2):end)));%second
    end
end
