% sandbox_SWR_Sub_oct2026

clear all


f_list = dir('pox*.mat');
k_idx = ones(length(f_list),1);

plot_flag =1;

% init
usr_n = {'ch','pos', 'sh' 'depth', 'fr', 'deep', 'rel_depth'};
% make a spike array
for ii = 1:length(usr_n)
    all_S.(usr_n{ii}) = [];
end
all_S.loc = [];
all_S.sess_idx = [];
all_S.sub_idx = [];
sess_id = {'TFCD1'  'TFCD2'  'TFCD3'  'TFCD4'  'TFCD5'};
sub_id = {'pox2217' 'pox3265'  'pox3567'  'pox3568'};

for iS = 1:length(f_list)
    % fname = 'pox3568_TFCD4'; % 3568 TFCD2/4/5 are all nice.


    % start with a nice session and generate some simple plots.

    this_sess = load([f_list(iS).name]);
    fname = fieldnames(this_sess.data);
    % Extract the session data for further analysis
    this_sess = this_sess.data.(fname{1}); % Access the first field of the session data

    % if length(strfind(fname{1}, '_')) > 1
    numText = regexp(fname{1}, '-?\d+\.?\d*', 'match', 'once');
    u_idx = strfind(fname{1}, '_');

    out.sess{iS} = upper(fname{1}(u_idx(end)+1:end));
    out.sub{iS} = ['pox' numText];

    if contains(out.sess{iS}, 'TL')
        out.sess{iS} = strrep(out.sess{iS}, 'TL', 'LT');
    end

    % usr_n = {'ch','pos', 'shank' 'depth', 'fr', 'deep', 'rel_depth'};
    % make a spike array
    for ii = 1:length(usr_n)
        all_S.(usr_n{ii}) = [all_S.(usr_n{ii}); this_sess.S.usr.(usr_n{ii})];
    end

    all_S.loc = [all_S.loc; this_sess.S.loc];
    % for jj = length(this_sess.S.t):-1:1
    %     isi(jj) = mode(diff(this_sess.S.t{jj}));
    % end

    all_S.sess_idx = [all_S.sess_idx; repmat(find(ismember(out.sess{iS}, sess_id)), length(this_sess.S.loc),1)];
    all_S.sub_idx = [all_S.sub_idx; repmat(find(ismember(out.sub{iS}, sub_id)), length(this_sess.S.loc),1)];

end

all_S.isi = (1./all_S.fr)*1000; 


% remove NaNs

nan_idx = isnan(all_S.fr); 
usr_n = fieldnames(all_S);
for ii = 1:length(usr_n)
    all_S.(usr_n{ii})(nan_idx) = []; 
end

pyr_idx = all_S.fr < 30; 

%% basic plots
c_ord = MS_linspecer(5);

figure(1)
clf

% subplot(2,2,[1 2 5 6])
% cla
% hold on
% scatter(all_S, all_S.isi)

subplot(2,2,1)
cla
hb = MS_bar_w_err(all_S.fr(all_S.loc==1 & pyr_idx), all_S.fr(all_S.loc==0 & pyr_idx),c_ord(1:2,:),1, 'ttest2', 1:2, 0); 
set(gca, 'YScale', 'linear','xticklabel', {'CA1' 'Sub'})
ylim([0 50])
ylabel('log Firing rate')
title ('Pyramidal Cells')

subplot(2,2,3)
cla
hb = MS_bar_w_err3(all_S.fr(all_S.loc==1 & all_S.deep==0 & pyr_idx), all_S.fr(all_S.loc==1 & all_S.deep==1 & pyr_idx), all_S.fr(all_S.loc==0 & pyr_idx),[c_ord(5,:); c_ord(3,:); c_ord(2,:)],1, 'anova1', 1:3); 
set(gca, 'YScale', 'linear','xticklabel', {'CA1_{deep}' 'CA1_{Super}' 'Sub'})
ylim([0 50])
ylabel('log Firing rate')

subplot(2,2,2)
cla
hb = MS_bar_w_err(all_S.fr(all_S.loc==1 & ~pyr_idx), all_S.fr(all_S.loc==0 & ~pyr_idx),c_ord(1:2,:),1, 'ttest2', 1:2, 0); 
set(gca, 'YScale', 'linear','xticklabel', {'CA1' 'Sub'})
ylim([0 100])
ylabel('log Firing rate')
title ('Fast Spiking Cells')

subplot(2,2,4)
cla
hb = MS_bar_w_err3(all_S.fr(all_S.loc==1 & all_S.deep==0 & ~pyr_idx), all_S.fr(all_S.loc==1 & all_S.deep==1 & ~pyr_idx), all_S.fr(all_S.loc==0 & ~pyr_idx),[c_ord(5,:); c_ord(3,:); c_ord(2,:)],1, 'anova1', 1:3); 
set(gca, 'YScale', 'linear','xticklabel', {'CA1_{deep}' 'CA1_{Super}' 'Sub'})
ylim([0 100])
ylabel('log Firing rate')


%% try some assembly stuff

