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


%%
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


this_rate = MS_spike2rate(this_sess.S, this_sess.csc.tvec, .025, 0); 


mov_ts = ts({this_sess.evts.t{find(ismember(this_sess.evts.label, 'mov'))}}, {'move_binary'});

mov_rate =  MS_spike2rate(mov_ts, this_rate.tvec, .025, 0);
move_idx =  mov_rate.data > 0; 

if length(move_idx) ~= length(this_rate.tvec)
    pad = length(move_idx) - length(this_rate.tvec);
    move_idx(end+(abs(pad))) = 0; 
end



this_asmbly = MS_asmbly_ephys(this_rate, move_idx);

% this_asmbly = MS_asmbly_members(this_asmbly)

%% get the phase locking to theta and maybe gamma

theta_csc = MS_select_tsd(this_sess.csc, 9);

mod_out = MS_phase_mod([], theta_csc, this_sess.S, move_idx, plot_flag)


%% check the phase locking of assemblies to theta




%% quick corr mat of the assemblies
corr_mat = NaN(size(this_asmbly.A_temp,2)); 

for ii = 1:size(this_asmbly.A_temp,2)
    for jj = 1:size(this_asmbly.A_temp,2)
        corr_mat(ii,jj) = corr(this_asmbly.A_proj(ii,:)', this_asmbly.A_proj(jj,:)'); 
    end
end

% replace diag with zero
% corr_mat = corr_mat - diag(diag(corr_mat)); 


d = 1 - corr_mat;
d(1:1+size(d,1):end) = 0;

Z = linkage(squareform(d), 'average');





figure(101)
% dendro
subplot(1,3,1)
[h, t, outperm] = dendrogram(Z);
title('Hierarchical Clustering of Correlation Matrix');
xlabel('Variable Index');
ylabel('Linkage Distance (1 - r)');

% sorted corr_mat
s_corr_mat = corr_mat(outperm,outperm); 
% replace lower mat with NaN; 
linear_idx = find(tril(ones(size(s_corr_mat))));
s_corr_mat(linear_idx) = NaN;


subplot(1,3,2)
imagesc(outperm, outperm,s_corr_mat)
xlabel('Assembly #')
ylabel('Assembly #')
c = colorbar;
c.Label.String = 'Corr';
clim([0 1]); 


subplot(1,3,3)
corr_mat(1:1+size(corr_mat,1):end) = 0;

cg = clustergram(corr_mat, 'Standardize', 'none', 'symmetric', 'false');

set(gca, 'YDir', 'reverse')



% try some clustering. 
% eva = evalclusters(corr_mat, 'kmeans', 'CalinskiHarabasz', 1:5)
%% plot it

figure(1010)
clf
this_csc = MS_select_tsd(this_sess.csc, {this_sess.swr_idx{1}});
MS_asmbly_ephys_raster(this_sess.S, this_rate.tvec, this_asmbly.A_temp, this_asmbly.A_proj,[1 5 6 7 8], this_csc, this_sess.S.loc)
