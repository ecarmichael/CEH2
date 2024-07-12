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
        
        S_dir{ii}(iD) = find(contains(S_loc,Days{iD}(end))); % Start location index
        D_ID{ii}(iD) = str2double(Days{iD}(end-1)); % day
        
    end
end



%% loop over sessions and get the corresponding position.

% maze specs

maze.cent = [958.3148 565.6295]; 
maze.rad = 423.4; 
maze.xv = maze.rad*1.1*cos(0:pi/100:2*pi)+maze.cent(1); 
maze.yv = maze.rad*1.1*sin(0:pi/100:2*pi)+maze.cent(2); 

c_ord = MS_linspecer(4);
lgd_names = [];

P_len = []; 
c_time = []; 
TQ = [];
RQ = [];
LQ = [];
OQ = []; 

sub = [3, 4, 5]; 
sess = 1:4; 
figure(101)
clf
for ii = 1:length(sub)
    t_day = [9 10 11 12]; 
    for iD = 1:length(t_day)
        %
        
        %         iD = 4;
        
        % find the label file for this subject/day.
        
        fname = dir(['labels*Day' num2str(D_ID{sub(ii)}(t_day(iD))) '*' M_ID{sub(ii)}(end-2:end) '*']);
        
        if isempty(fname); continue; end
        
        pos = MS_SLEAP2TSD(fname.name, 30);
        
        pos_r = restrict(pos, pos.tvec(tstart{sub(ii)}(t_day(iD))),  pos.tvec(tend{sub(ii)}(t_day(iD))));
        
%         % fill in missing values outside of the maze. 
%         d_tbl = reshape(1:size(pos_r.data,1), 2, size(pos_r.data,1)/2);
%         for iF = 1:length(pos_r.label)-2
%             
%             [isin, ison] = inpolygon(pos_r.data(d_tbl(iF),:), pos_r.data(d_tbl(iF+1),:), maze.xv, maze.yv);
%             pos_r.data(d_tbl(iF),~(isin|ison)) = NaN; 
%             pos_r.data(d_tbl(iF+1),~(isin|ison)) = NaN; 
%             
%             pos_r.data(d_tbl(iF),:) = fillmissing(pos_r.data(d_tbl(iF)), 'spline', 'EndValues','nearest');
%             pos_r.data(d_tbl(iF+1),:) = fillmissing(pos_r.data(d_tbl(iF+1)), 'spline', 'EndValues','nearest');
% 
%         end
        
        %%
        subplot(2, length(sub), 1:length(sub));
        hold on
        if ii ==1 && iD == 1
            plot(pos.tvec, pos.data([1 3 5 7],:), 'k')
        end
        plot(pos_r.tvec(1), pos_r.data([1 3 5 7],1), 'd','color', c_ord(iD,:), 'MarkerSize', 20)
        plot(pos_r.tvec, pos_r.data([1 3 5 7],:),'color', c_ord(iD,:))
        plot(pos_r.tvec(end), pos_r.data([1 3 5 7],end), '*','color', c_ord(iD,:), 'MarkerSize', 20)
        
        % plot check
        subplot(2, length(sub), ii + length(sub));
        hold on
        if iD ==1
            imagesc(pos.mean_frame); colormap('gray')
            plot(maze.xv, maze.yv)
        end
        plot(pos_r.data(1,:), pos_r.data(2,:),'.', 'color', c_ord(iD,:))
        plot(pos_r.data(3,:), pos_r.data(4,:),'.', 'color', c_ord(iD,:)*.8,  'MarkerSize', 4)
        
        plot(nanmean(pos_r.data(1,1:10)), nanmean(pos_r.data(2,1:10)),'d', 'color', c_ord(iD,:), 'MarkerSize', 20)
        fprintf('Trial: <strong>%s</strong> nan frames %0.0f /%0.0f (%0.2f%%)\n', S_loc{S_dir{sub(ii)}(iD)}, sum(isnan(pos_r.data(5,:))), length(pos_r.data(5,:)), sum(isnan(pos_r.data(5,:)))/length(pos_r.data(5,:)))
        lgd_names{end+1} =  S_loc{S_dir{sub(ii)}(t_day(iD))};
        lgd_names{end+1} =  [S_loc{S_dir{sub(ii)}(t_day(iD))} 'start'];
        title(M_ID{sub(ii)})
        
        %% compute simple measures
        
        P_len{sub(ii)}(t_day(iD)) = nansum(sqrt(diff(pos_r.data(1,:)).^2 + diff(pos_r.data(2,:)).^2)); 
        c_time{sub(ii)}(t_day(iD)) = pos_r.tvec(end) - pos_r.tvec(1); 
        
        % divide the maze into quadrants. 
        
        QT_idx =(pos_r.data(1,:) > maze.cent(1)) & (pos_r.data(2,:) < maze.cent(2));
        QO_idx =(pos_r.data(1,:) < maze.cent(1)) & (pos_r.data(2,:) > maze.cent(2));
        QR_idx =(pos_r.data(1,:) < maze.cent(1)) & (pos_r.data(2,:) < maze.cent(2));
        QL_idx =(pos_r.data(1,:) > maze.cent(1)) & (pos_r.data(2,:) > maze.cent(2));
        
        % % time in each 
        QT{sub(ii)}(t_day(iD)) = sum(QT_idx/length(QT_idx))*100; 
        QO{sub(ii)}(t_day(iD)) = sum(QO_idx/length(QO_idx))*100;
        QR{sub(ii)}(t_day(iD)) = sum(QR_idx/length(QR_idx))*100; 
        QL{sub(ii)}(t_day(iD)) = sum(QL_idx/length(QL_idx))*100; 


        figure(191)
        hold on
        plot(pos_r.data(1,QT_idx), pos_r.data(2, QT_idx), '.r'); 
        plot(pos_r.data(1,QO_idx), pos_r.data(2, QO_idx), '.b'); 
        plot(pos_r.data(1,QR_idx), pos_r.data(2, QR_idx), '.g'); 
        plot(pos_r.data(1,QL_idx), pos_r.data(2, QL_idx), '.y'); 

        
        
    end
end

legend(lgd_names)

% end



% format as a table

Data_tbl = [];
% Data_tbl.M_ID =

