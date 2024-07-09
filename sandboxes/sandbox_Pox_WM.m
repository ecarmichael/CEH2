% % sandbox_WM_split


cd('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\PoxR1\WM'); 


POX_tbl = readtable("POX_WM - Sheet1.csv");

M_idx = find(contains(POX_tbl.ID, 'Pox')); 
Days = POX_tbl.Properties.VariableNames(2:end); 

S_loc = {'N', 'E', 'S', 'W'};
Fs = 30; 


for ii = length(M_idx):-1:1
    M_ID{ii} = POX_tbl.ID{M_idx(ii)};
    
    for iD = length(Days):-1:1
        time{ii}(iD) = POX_tbl.(Days{iD})(M_idx(ii)); % get the watch time. 
        tstart{ii}(iD) = POX_tbl.(Days{iD})(M_idx(ii)+1); % get the start frame
        
        if time{ii}(iD) == 60  % if it was a timeout then set to 60s
           tend{ii}(iD) =  tstart{ii}(iD)+(60*Fs); 
        else
        tend{ii}(iD) = POX_tbl.(Days{iD})(M_idx(ii)+2); % get the end frame. 
        end
        
        f_time{ii}(iD) = (tend{ii}(iD) - tstart{ii}(iD))./Fs; % frame clock time. 
        
        S_dir{ii}(iD) = find(contains(Days{iD}(end), S_loc)); % Start location index
        D_ID{ii}(iD) = str2double(Days{iD}(end-1)); % day 
        
    end
end

% format as a table

Data_tbl = []; 
% Data_tbl.M_ID = 

