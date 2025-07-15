% % sandbox_WM_split

if ispc
    cd('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\PoxR1\WM');
elseif ismac
    cd('/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/WM');
else
    cd('/home/swlab/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/WM'); 
end




POX_tbl = readtable("POX_WM - Sheet1.csv");

M_idx = find(contains(POX_tbl.ID, 'Pox'));
Days = POX_tbl.Properties.VariableNames(2:end-2);

S_loc = {'N', 'E', 'S', 'W'};
P_loc = {'NE', 'SE', 'SW', 'NW'};
Fs = 30;


for ii = length(M_idx):-1:1
    M_ID{ii} = POX_tbl.ID{M_idx(ii)};
    
    for iD = length(Days):-1:1
        u_time{ii}(iD) = POX_tbl.(Days{iD})(M_idx(ii)); % get the watch time.
        tstart{ii}(iD) = POX_tbl.(Days{iD})(M_idx(ii)+1); % get the start frame
        
        if u_time{ii}(iD) == 60  % if it was a timeout then set to 60s
            tend{ii}(iD) =  tstart{ii}(iD)+(60*Fs);
        else
            tend{ii}(iD) = POX_tbl.(Days{iD})(M_idx(ii)+2); % get the end frame.
        end
        
        f_time{ii}(iD) = (tend{ii}(iD) - tstart{ii}(iD))./Fs; % frame clock time.
        
        if contains(lower(Days{iD}(1:5)), lower('Probe'))
            S_dir{ii}(iD) =find(contains(P_loc,Days{iD}(end-1:end))); % Start location index
            D_ID{ii}(iD) = 999; % day
        else

            S_dir{ii}(iD) = find(contains(S_loc,Days{iD}(end))); % Start location index
            D_ID{ii}(iD) = str2double(Days{iD}(end-1)); % day
        end
    end
end



%% loop over sessions and get the corresponding position.

all_data = [];

% maze specs
maze.cent = [958.3148 565.6295];
maze.rad = 423.4;
maze.xv = maze.rad*1.1*cos(0:pi/100:2*pi)+maze.cent(1);
maze.yv = maze.rad*1.1*sin(0:pi/100:2*pi)+maze.cent(2);
maze.platform_diameter = 11.5; 
maze.platform_cent = [1134.3 367.68];
maze.units = 'cm';
maze.diameter_cm = 115; 
maze.diameter_px = 870; 
maze.conv_fact = maze.diameter_px/maze.diameter_cm; 

c_ord = MS_linspecer(4);
lgd_names = [];

P_len = [];
c_time = [];
TQ = [];
RQ = [];
LQ = [];
OQ = [];

sub = 1:length(M_ID);
sess = 1:4;

figure(101)
clf
for ii = 1:length(sub)
    
    %%
    t_day = 1:length(D_ID{sub(ii)});
    for iD = 1:length(t_day)
        %
        
        %         iD = 4;
        
        % find the label file for this subject/day.
        
        fname = dir(['labels*Day' num2str(D_ID{sub(ii)}(t_day(iD))) '*' M_ID{sub(ii)}(end-2:end) '*']);
        
        if isempty(fname); continue; end
        
        pos = MS_SLEAP2TSD(fname.name, 30, [], [maze.conv_fact maze.conv_fact]);
        
        pos_r = restrict(pos, pos.tvec(tstart{sub(ii)}(t_day(iD))),  pos.tvec(tend{sub(ii)}(t_day(iD))));
        
        % center the data to the middle of the pool;
        for jj = 1:2:size(pos_r.data,1)-2
            pos_r.data(jj,:) = pos_r.data(jj,:) - maze.cent(1)/maze.conv_fact; 
            pos_r.data(jj+1,:) = pos_r.data(jj+1,:) - maze.cent(2)/maze.conv_fact; 
        end
        
        %%
        figure(191)
        
        subplot(2, length(sub), 1:length(sub));
        hold on
        if ii ==1 && S_dir{sub(ii)}(t_day(iD)) == 1
            plot(pos.tvec, pos.data([1 3 5 7],:), 'k')
        end
        plot(pos_r.tvec(1), pos_r.data([1 3 5 7],1), 'd','color', c_ord(S_dir{sub(ii)}(t_day(iD)),:), 'MarkerSize', 20)
        plot(pos_r.tvec, pos_r.data([1 3 5 7],:),'color', c_ord(S_dir{sub(ii)}(t_day(iD)),:))
        plot(pos_r.tvec(end), pos_r.data([1 3 5 7],end), '*','color', c_ord(S_dir{sub(ii)}(t_day(iD)),:), 'MarkerSize', 20)
        
        % plot check
        subplot(2, length(sub), ii + length(sub));
        hold on
        set(gca, 'ydir', 'reverse')
        if iD ==1
            x = 1:size(pos.mean_frame,2); 
            y = 1:size(pos.mean_frame,1);
            x  = x./maze.conv_fact; 
            y  = y./maze.conv_fact;

            x = x - (maze.cent(1)/maze.conv_fact); 
            y = y - (maze.cent(2)/maze.conv_fact); 


            imagesc(x, y,pos.mean_frame); colormap('gray')
            plot((maze.xv./maze.conv_fact)- (maze.cent(1)/maze.conv_fact), (maze.yv./maze.conv_fact)- (maze.cent(2)/maze.conv_fact), '.w')
        end
        plot(pos_r.data(1,:), pos_r.data(2,:),'.', 'color', c_ord(S_dir{sub(ii)}(t_day(iD)),:))
        plot(pos_r.data(3,:), pos_r.data(4,:),'.', 'color', c_ord(S_dir{sub(ii)}(t_day(iD)),:)*.8,  'MarkerSize', 4)
        
        plot(nanmean(pos_r.data(1,1:10)), nanmean(pos_r.data(2,1:10)),'d', 'color', c_ord(S_dir{sub(ii)}(t_day(iD)),:), 'MarkerSize', 20)
        fprintf('Trial: <strong>%s</strong> nan frames %0.0f /%0.0f (%0.2f%%)\n', S_loc{S_dir{sub(ii)}(iD)}, sum(isnan(pos_r.data(5,:))), length(pos_r.data(5,:)), sum(isnan(pos_r.data(5,:)))/length(pos_r.data(5,:)))
        lgd_names{end+1} =  S_loc{S_dir{sub(ii)}(t_day(iD))};
        lgd_names{end+1} =  [S_loc{S_dir{sub(ii)}(t_day(iD))} 'start'];
        title(M_ID{sub(ii)})
        
        %% compute simple measures
        
        P_len{sub(ii)}(t_day(iD)) = nansum(sqrt(diff(pos_r.data(1,:)).^2 + diff(pos_r.data(2,:)).^2));
        c_time{sub(ii)}(t_day(iD)) = pos_r.tvec(end) - pos_r.tvec(1);
        
        % divide the maze into quadrants.
        
        QT_idx =(nanmean(pos_r.data([1 3 5],:)) > maze.cent(1)) & (nanmean(pos_r.data([2 4 6],:)) < maze.cent(2));
        QO_idx =(nanmean(pos_r.data([1 3 5],:)) < maze.cent(1)) & (nanmean(pos_r.data([2 4 6],:)) > maze.cent(2));
        QR_idx =(nanmean(pos_r.data([1 3 5],:)) < maze.cent(1)) & (nanmean(pos_r.data([2 4 6],:)) < maze.cent(2));
        QL_idx =(nanmean(pos_r.data([1 3 5],:)) > maze.cent(1)) & (nanmean(pos_r.data([2 4 6],:)) > maze.cent(2));
        
        % % time in each
        QT{sub(ii)}(t_day(iD)) = sum(QT_idx/length(QT_idx))*100;
        QO{sub(ii)}(t_day(iD)) = sum(QO_idx/length(QO_idx))*100;
        QR{sub(ii)}(t_day(iD)) = sum(QR_idx/length(QR_idx))*100;
        QL{sub(ii)}(t_day(iD)) = sum(QL_idx/length(QL_idx))*100;
        
        
        figure(191)
        hold on
        plot(nanmean(pos_r.data([1 3 5],QT_idx)), nanmean(pos_r.data([2  4 6], QT_idx)), '.r');
        plot(nanmean(pos_r.data([1 3 5],QO_idx)), nanmean(pos_r.data([2  4 6], QO_idx)), '.b');
        plot(nanmean(pos_r.data([1 3 5],QR_idx)), nanmean(pos_r.data([2  4 6], QR_idx)), '.g');
        plot(nanmean(pos_r.data([1 3 5],QL_idx)), nanmean(pos_r.data([2  4 6], QL_idx)), '.y');
        
        % save the data for later.
        all_data.(M_ID{sub(ii)}).(Days{iD}).maze = maze;
        all_data.(M_ID{sub(ii)}).(Days{iD}).pos_r = pos_r;
        
        
    end
end

% legend(lgd_names)

% end

%% PRobe only
n = ceil(length(M_ID)/4)+1;
m = 3;


for ii = 1:length(M_ID)
    for iD = length(D_ID{ii})
        
        
        fname = dir(['labels*Probe*' M_ID{ii}(end-2:end) '*']);
        
        if isempty(fname); continue; end
        
        pos = MS_SLEAP2TSD(fname.name, 30);
        
        pos_r = restrict(pos, pos.tvec(tstart{ii}(iD)),  pos.tvec(tend{ii}(iD)));
        
        
        QT_idx =(nanmean(pos_r.data([1 3 5],:)) > maze.cent(1)) & (nanmean(pos_r.data([2 4 6],:)) < maze.cent(2));
        QO_idx =(nanmean(pos_r.data([1 3 5],:)) < maze.cent(1)) & (nanmean(pos_r.data([2 4 6],:)) > maze.cent(2));
        QR_idx =(nanmean(pos_r.data([1 3 5],:)) < maze.cent(1)) & (nanmean(pos_r.data([2 4 6],:)) < maze.cent(2));
        QL_idx =(nanmean(pos_r.data([1 3 5],:)) > maze.cent(1)) & (nanmean(pos_r.data([2 4 6],:)) > maze.cent(2));
        
        % % time in each
        QT{ii}(iD) = sum(QT_idx/length(QT_idx))*100;
        QO{ii}(iD) = sum(QO_idx/length(QO_idx))*100;
        QR{ii}(iD) = sum(QR_idx/length(QR_idx))*100;
        QL{ii}(iD) = sum(QL_idx/length(QL_idx))*100;
        
        
        subplot(m, n, ii)
        hold on
        set(gca, 'ydir', 'reverse')
        
        imagesc(pos.mean_frame);
        colormap('parula')
        scatter(nanmean(pos_r.data([1, 3, 5],:)), nanmean(pos_r.data([2, 4, 6],:)),5, 1:length(pos_r.tvec), '.')
        title({strrep(M_ID{ii}, '_', ' ') ; ['T: ' num2str(QT{ii}(iD),3) ' |O: ' num2str(QO{ii}(iD),3) ' |O: ' num2str(QR{ii}(iD),3) ' |L: ' num2str(QL{ii}(iD),3)]})
        axis square
        axis off
        P_QT(ii) = QT{ii}(iD);
        P_QO(ii) = QO{ii}(iD);
        P_QR(ii) = QR{ii}(iD);
        P_QL(ii) = QL{ii}(iD);
    end
end

% MS_bar(mean([P_QT, P_QO, )

%% Convert to 'watermaze.csv' for PATHFINDER
subs = fieldnames(all_data);

target = {'Nose','Head'};


for iS =  1:length(subs)
    %%
    
%     figure(iS)
%     clf
    days= fieldnames(all_data.(subs{iS}));
    tbl = table();
    
    % get the max number of data points for the table
    max_d = [];
    for iD = length(days):-1:1
        max_d(iD) = length(all_data.(subs{iS}).(days{iD}).pos_r.tvec);
    end
    d_ord = MS_linspecer(length(days));
    
    for iD = 1:length(days)
        
        d_idx = strfind(all_data.(subs{iS}).(days{iD}).pos_r.cfg.fname{1}, '2024');
        this_date = all_data.(subs{iS}).(days{iD}).pos_r.cfg.fname{1}(d_idx:d_idx+9);
        this_date = datetime(this_date, 'InputFormat', 'yyyy_MM_dd', 'Format', 'yyyy/MM/dd');
        
        %         this_x = fillmissing(mean(all_data.(subs{iS}).(days{iD}).pos_r.data([1, 3, 5 ],:))', 'spline');
        %         this_y = fillmissing(mean(all_data.(subs{iS}).(days{iD}).pos_r.data([2, 4, 6],:))', 'spline');
        this_x = mean(all_data.(subs{iS}).(days{iD}).pos_r.data([1, 3, 5 ],:))' - all_data.(subs{iS}).(days{iD}).maze.cent(1);
        this_y = mean(all_data.(subs{iS}).(days{iD}).pos_r.data([2, 4, 6],:))'- all_data.(subs{iS}).(days{iD}).maze.cent(2);
        this_t = all_data.(subs{iS}).(days{iD}).pos_r.tvec';
        
        % find the first non-nan and remove the data up to the first data
        % point.
        
        dat_idx = find(~isnan(this_x));
        
        if dat_idx(1) >1
            this_x(1:dat_idx(1)) = [];
            this_y(1:dat_idx(1)) = [];
            this_t(1:dat_idx(1)) = [];
        end
        
        this_t = this_t - this_t(1);
        %         Flip the Y so that the platform is in the top right quadrant.
        t_y = this_y;
        t_y(this_y>0) = -this_y(this_y>0);
        t_y(this_y<0) = abs(this_y(this_y<0));
        
        this_y = t_y;
        
        if length(this_x) <= max(max_d)
            this_x(length(this_x)+1:max(max_d)) = NaN;
            this_y(length(this_y)+1:max(max_d)) = NaN;
            this_t(length(this_t)+1:max(max_d)) = NaN;
        end
        
        
        tbl = addvars(tbl, [days{iD}; 'x'; num2cell(this_x)]);
        tbl = addvars(tbl, [char(this_date); 'y'; num2cell(this_y)]);
        tbl = addvars(tbl, [' '; 't'; num2cell(this_t)]);
        
        hold on
        plot(this_x, this_y, '.', 'color', d_ord(iD,:));
    end
    
    % try replacing NaNs by converting to a cell arraty
    c_tbl = table2cell(tbl);
    
    % mask = cellfun(@(C), )
    
    for ii = size(c_tbl,1):-1:1
        nan_idx = cell2mat(cellfun(@isnan, c_tbl(ii,:), 'UniformOutput',false));
        
        % c_tbl(ii,3:end) = num2str(c_tbl(ii,3:end));
        c_tbl(ii,nan_idx) = {' '};
    end
    
    tbl_out = cell2table(c_tbl);
    
    %write a .csv per animal.
    writetable(tbl, [subs{iS} '_WM.csv'], 'WriteVariableNames',1)
    % tbl_in = readtable([subs{iS} '_WM.csv'])
    
end








