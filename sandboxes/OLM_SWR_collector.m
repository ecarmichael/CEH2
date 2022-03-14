%% OLM_SWR_collector

%% process all sessions
data_dir = '/home/williamslab/Dropbox (Williams Lab)/Data SWR sm';

cd(data_dir)
sess_list = dir('*day*');

for ii = 1:length(sess_list)
   cd(sess_list(ii).name)
   
   OLM_SWR_screener(30); 
   close all
   cd(data_dir)

end

%% get the data

data = []; all_label = [];
data_t = [];
inter_dir = '/home/williamslab/Dropbox (Williams Lab)/Data SWR sm/inter';
cd(inter_dir)
sessions = dir('*.mat'); 


for iS = 1:length(sessions)
    
   
    load([inter_dir filesep sessions(iS).name]); 
    
    all_data{iS}.sub = data_out.f_info.subject;
    all_data{iS}.sess = data_out.f_info.session;

    all_SWRs{iS} = data_out.SWR_evts; 
    
    
    % get the SWR durations per time bin;
    evt_dur = [];
    for ii  = length(data_out.t):-1:1
        if ii == length(data_out.t)
                    this_swr = restrict(data_out.SWR_evts, data_out.t(ii), data_out.csc.tvec(end)); 
        else
                    this_swr = restrict(data_out.SWR_evts, data_out.t(ii), data_out.t(ii+1)); 
        end
        evt_dur(ii) = mean(this_swr.tend - this_swr.tstart);
    end

%     all_data{iS}.sub = data_out.
    
    if isfield(data_out, 't_minus') 
        all_data{iS}.t = [((data_out.t_minus - data_out.t_minus(1))/60/60)-1 data_out.t_zero]; 
        all_data{iS}.rate = [data_out.pre_nEvts./ data_out.pre_nDur,  data_out.nEvts./ data_out.nDur];
        data_len(iS,:) = evt_dur;
        data(iS,:) = all_data{iS}.rate(1:11); 
        data_t(iS,:) = [(data_out.nEvts./ data_out.nDur) ./ nanmean(data_out.pre_nEvts./ data_out.pre_nDur)]; 
    else
        all_data{iS}.t = [nan(1,3), data_out.t_zero]; 
        all_data{iS}.rate = [nan(1,3),  data_out.nEvts./ data_out.nDur];
        data_len(iS,:) =  evt_dur; 

        data(iS,:) = all_data{iS}.rate(1:11); 
    end

    all_label(iS,2) = str2double(data_out.f_info.session(end));
    if strcmpi(data_out.f_info.subject, 'JB1556')
        all_label(iS,1) = 1;
    elseif strcmpi(data_out.f_info.subject, 'JB1446')
        all_label(iS,1) = 2;
    end
    
    data_t(data_t == 0) = nan; 
    
end
%% plot

% this_data = data_t; 

tvec = all_data{1}.t(3:end)+(mode(diff(all_data{1}.t(4:end)))/2);
tvec = tvec(1:end-1)

this_data = data_len(:,1:end-1);
% tvec = all_data{1}.t+(mode(diff(all_data{1}.t))/2);
c_ord = linspecer(3); 
cp = [1,3, 5]; 
subj = [1,2]; 

cp_ix = ismember(all_label(:,2),  cp); 
s1_cp = (all_label(:,1) == 1) & cp_ix; 
s2_cp = (all_label(:,1) == 2) & ~cp_ix; 
s1_v = (all_label(:,1) == 1) & ~cp_ix; 
s2_v = (all_label(:,1) == 2) & cp_ix;

figure(301)
clf
hold on
% subj 1 cp21
h1 = shadedErrorBar(tvec, nanmean(this_data(s1_cp,:), 1), nanstd(this_data(s1_cp,:),[],1) / sqrt(length(this_data(s1_cp,:))));
h1.mainLine.LineWidth = 3; h1.mainLine.Color = c_ord(1,:); 

h1.edge(1).LineWidth = 1; h1.edge(1).Color = [c_ord(1,:) , .3]; 
h1.edge(2).LineWidth = 1; h1.edge(2).Color = [c_ord(1,:) , .3]; 
 
h1.patch.FaceColor = c_ord(1,:); h1.patch.FaceAlpha = .2;

% subj 2 cp21
h2 = shadedErrorBar(tvec, nanmean(this_data(s2_cp,:), 1), nanstd(this_data(s2_cp,:),[],1) / sqrt(length(this_data(s2_cp,:))));
h2.mainLine.LineWidth = 3; h2.mainLine.Color = c_ord(2,:); 

h2.edge(1).LineWidth = 1; h2.edge(1).Color = [c_ord(2,:) , .3]; 
h2.edge(2).LineWidth = 1; h2.edge(2).Color = [c_ord(2,:) , .3]; 
 
h2.patch.FaceColor = c_ord(2,:); h2.patch.FaceAlpha = .2;

% add in the vehicle days
h3 = shadedErrorBar(tvec, nanmean(this_data(s1_v,:), 1), nanstd(this_data(s1_v,:),[],1) / sqrt(length(this_data(s1_v,:))));
h3.mainLine.LineWidth = 3; h3.mainLine.Color = c_ord(1,:); h3.mainLine.LineStyle = '--'; 

h3.edge(1).LineWidth = 1; h3.edge(1).Color = [c_ord(1,:) , .3]; h3.edge(1).LineStyle = '--'; 
h3.edge(2).LineWidth = 1; h3.edge(2).Color = [c_ord(1,:) , .3]; h3.edge(2).LineStyle = '--'; 
 
h3.patch.FaceColor = c_ord(1,:); h3.patch.FaceAlpha = .2;

% subj 2 cp21
h4 = shadedErrorBar(tvec, nanmean(this_data(s2_v,:), 1), nanstd(this_data(s2_v,:),[],1) / sqrt(length(this_data(s2_v,:))));
h4.mainLine.LineWidth = 3; h4.mainLine.Color = c_ord(2,:); h4.mainLine.LineStyle = '--'; 

h4.edge(1).LineWidth = 1; h4.edge(1).Color = [c_ord(2,:) , .3]; h4.edge(1).LineStyle = '--'; 
h4.edge(2).LineWidth = 1; h4.edge(2).Color = [c_ord(2,:) , .3]; h4.edge(2).LineStyle = '--'; 
 
h4.patch.FaceColor = c_ord(2,:); h4.patch.FaceAlpha = .2;

h = get(gca, 'Children');
legend(h(fliplr([1 5 9 13])),{'1556-cp21','1446-cp21', '1556-v','1446-v'} )
% legend({'', '', '','1556-cp21','', '', '', '1446-cp21','', '', '', '1556-v','', '', '', '1446-v'})


%% bar plot version
figure(401)

hold on;
% h = bar(tvec, [nanmean(this_data(s1_cp,:));nanmean(this_data(s1_v,:))]',  'FaceColor', 'flat', 'LineWidth', .9);
% h(1).CData(1, :) = cl(1, :);
% h(1).CData(2, :) = cl(2, :);

% e = errorbar(1:length([nanmean(this_data(s1_cp,:));nanmean(this_data(s2_cp,:))]'), [nanmean(this_data(s1_cp,:));nanmean(this_data(s2_cp,:))], [nanstd(this_data(s2_cp,:),[],1) / sqrt(length(this_data(s2_cp,:))) ; nanstd(this_data(s2_cp,:),[],1) / sqrt(length(this_data(s2_cp,:)))']', '.k', 'LineWidth', .9);
% errorbar_groups(nanmean(this_data(s1_cp,:), 1), nanstd(this_data(s1_cp,:),[],1) / sqrt(length(this_data(s1_cp,:))))

boxch


%% raincloud

% format the data in to data, measure,group

% time
for ii = 1:length(this_data)
    
    for jj = 1:size(this_data,1)
        
        
        
    end
end

rm_raincloud(this_data(s1_cp,:)')

