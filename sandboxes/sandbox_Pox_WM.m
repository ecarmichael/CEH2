% % sandbox_WM_split

if ispc
cd('C:\Users\ecarm\Williams Lab Dropbox\Williams Lab Team Folder\Eric\PoxR1\WM');

elseif ismac
    cd('/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/PoxR1/WM');
end




POX_tbl = readtable("POX_WM - Sheet1.csv");

M_idx = find(contains(POX_tbl.ID, 'Pox'));
Days = POX_tbl.Properties.VariableNames(2:end-2);

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
c_ord = MS_linspecer(4);
lgd_names = [];

sub = [3, 4, 5]; 
sess = 1:4; 
figure(101)
clf
for ii = 1:length(sub)
    %
    for iD = 1:
        %
        
        %         iD = 4;
        
        % find the label file for this subject/day.
        
        fname = dir(['labels*Day' num2str(D_ID{sub(ii)}(iD)) '*' M_ID{sub(ii)}(end-2:end) '*']);
        
        if isempty(fname); continue; end
        
        pos = MS_SLEAP2TSD(fname.name, 30);
        
        pos_r = restrict(pos, pos.tvec(tstart{sub(ii)}(iD)),  pos.tvec(tend{sub(ii)}(iD)));
        
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
        end
        plot(pos_r.data(1,:), pos_r.data(2,:),'.', 'color', c_ord(iD,:))
        plot(pos_r.data(3,:), pos_r.data(4,:),'.', 'color', c_ord(iD,:)*.2)

        plot(nanmean(pos_r.data(1,1:10)), nanmean(pos_r.data(2,1:10)),'d', 'color', c_ord(iD,:), 'MarkerSize', 20)
        fprintf('Trial: <strong>%s</strong> nan frames %0.0f /%0.0f (%0.2f%%)\n', S_loc{S_dir{sub(ii)}(iD)}, sum(isnan(pos_r.data(5,:))), length(pos_r.data(5,:)), sum(isnan(pos_r.data(5,:)))/length(pos_r.data(5,:)))
        lgd_names{end+1} =  S_loc{S_dir{sub(ii)}(iD)};
        lgd_names{end+1} =  [S_loc{S_dir{sub(ii)}(iD)} 'start'];
        title(M_ID{sub(ii)})
        
        
        
    end
end

legend(lgd_names)

% end



% format as a table

Data_tbl = [];
% Data_tbl.M_ID =

