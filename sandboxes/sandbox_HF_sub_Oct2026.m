%% sandbox_HF_Process_SWR


if ismac
    data_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Pox/';
    inter_dir = '/Users/ecar/Williams Lab Dropbox/Williams Lab Team Folder/Eric/Wheel/Inter_data/HF_Sub_SWR/inter_offset';
elseif ispc
    data_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\Pox\';
    inter_dir = 'C:\Users\ecar\Williams Lab Dropbox\Williams Lab Team Folder\Eric\Wheel\Inter_data\HF_Sub_SWR\inter_offset';
end

fname = dir([inter_dir filesep '*.mat']);

for iS = 1:length(fname)

    load([fname(iS).folder filesep fname(iS).name])

    this_name = fieldnames(data);

    this_sess = data.(this_name{1});
    parts = strsplit(this_name{1}, '_');
    sub_name = parts{1};
    sess_name = parts{end};

    % sess_idx = contains(all_csc_dir, sess_name);
    % sub_idx = contains(lower(all_csc_dir), sub_name);

    % this_idx = find(sess_idx & sub_idx);

    % fprintf('Inter: %s   =   %s   \n', this_name{1}, all_csc_dir{this_idx})



    if isfield(this_sess, 'swrs_ca1')
        disp(length(this_sess.swrs_ca1.tstart))


      % session
        sub_count{iS} = sub_name; 
        sess_count{iS} = sess_name; 

        % spikes
        n_Ca1_py(iS) = sum(this_sess.S.loc & this_sess.S.usr.fr <20); 
        n_Ca1_int(iS) = sum(this_sess.S.loc & this_sess.S.usr.fr >=20); 

        n_sub_py(iS) = sum(this_sess.S.loc==0 & this_sess.S.usr.fr <20); 
        n_sub_int(iS) = sum(this_sess.S.loc==0 & this_sess.S.usr.fr >=20); 


        %SWR
        swr_count(iS,1) = length(this_sess.swrs_ca1.tstart); 
        swr_count(iS,2) = length(this_sess.swrs_sub.tstart);

        % pre rate
        csc_r = restrict(this_sess.csc, this_sess.mov_iv); 
        csc_r_p = restrict(this_sess.csc, this_sess.csc.tvec(1), this_sess.csc.tvec(1)+300); 
        pre_ca1_swrs = restrict(this_sess.swrs_ca1, this_sess.csc.tvec(1), this_sess.csc.tvec(1)+300); 
        swr_ca1_rate(iS) = length(pre_ca1_swrs.tstart) / (length(csc_r_p.tvec)/csc_r_p.cfg.hdr{1}.SamplingFrequency); 

        pre_sub_swrs = restrict(this_sess.swrs_sub, this_sess.csc.tvec(1), this_sess.csc.tvec(1)+300); 
        swr_sub_rate(iS) = length(pre_sub_swrs.tstart)  / (length(csc_r_p.tvec)/csc_r_p.cfg.hdr{1}.SamplingFrequency); 


        % post rate
            csc_r_p = restrict(this_sess.csc, this_sess.csc.tvec(end) - 300, this_sess.csc.tvec(end)); 

        post_ca1_swrs = restrict(this_sess.swrs_ca1, this_sess.csc.tvec(end) - 300, this_sess.csc.tvec(end));
        swr_ca1_rate(iS) = length(post_ca1_swrs.tstart)  / (length(csc_r_p.tvec)/csc_r_p.cfg.hdr{1}.SamplingFrequency); 

        post_sub_swrs = restrict(this_sess.swrs_sub, this_sess.csc.tvec(end) - 300, this_sess.csc.tvec(end));
        swr_sub_rate(iS) = length(post_sub_swrs.tstart)  / (length(csc_r_p.tvec)/csc_r_p.cfg.hdr{1}.SamplingFrequency); 

        %

    end
end

keep_idx = swr_count(:,1) > 100; 
keep_idx(2) = 0; 
%% quick plot of the number per session x subject

c_ord = MS_linspecer(5); 

sub_n = unique(sub_count); 

sub_id = zeros(size(sub_count));
for ii = 1:length(sub_n)
    sub_id(ismember(sub_count, sub_n{ii})) = ii;
end

figure(101)
clf

subplot(2,2,1)
scatter(sub_id(keep_idx)+MS_randn_range(1,length(sub_id(keep_idx)),-.1, .1), swr_count(keep_idx,1), 50, 'markerfacecolor', c_ord(1,:))
ylim([0 max(swr_count(keep_idx,:), [], 'all')])
ylabel('basline SWR count')
set(gca,'xtick',1:length(sub_n),  'XTickLabel', sub_n)

subplot(2,2,2)
scatter(sub_id(keep_idx)+MS_randn_range(1,length(sub_id(keep_idx)),-.1, .1), swr_count(keep_idx,2), 50, 'markerfacecolor', c_ord(2,:))
set(gca, 'XTickLabel', sub_n)
ylim([0 max(swr_count, [], 'all')])
set(gca,'xtick',1:length(sub_n),  'XTickLabel', sub_n)


subplot(2,2,3)
scatter(sub_id(keep_idx)+MS_randn_range(1,length(sub_id(keep_idx)),-.1, .1), swr_ca1_rate(keep_idx), 50, 'markerfacecolor', c_ord(1,:))
ylim([0 max([swr_ca1_rate(keep_idx) swr_sub_rate(keep_idx)])])
ylabel('basline SWR rate (Hz)')
set(gca,'xtick',1:length(sub_n),  'XTickLabel', sub_n)

subplot(2,2,4)
scatter(sub_id(keep_idx)+MS_randn_range(1,length(sub_id(keep_idx)),-.1, .1), swr_sub_rate(keep_idx), 50, 'markerfacecolor', c_ord(2,:))
set(gca, 'XTickLabel', sub_n)
ylim([0 max([swr_ca1_rate(keep_idx) swr_sub_rate(keep_idx)])])
set(gca,'xtick',1:length(sub_n),  'XTickLabel', sub_n)



figure(102)
clf
subplot(2,2,1)
hold on
boxchart(sub_id(keep_idx), n_Ca1_py(keep_idx), 'BoxFaceColor',c_ord(1,:)); 
scatter(sub_id(keep_idx)+MS_randn_range(1,length(sub_id(keep_idx)),-.1, .1), n_Ca1_py(keep_idx), 50, 'markerfacecolor', c_ord(1,:), 'markeredgecolor', c_ord(1,:))
ylim([0 max([n_Ca1_py(keep_idx) n_sub_py(keep_idx)])])
ylabel('N Pyramidal Cells')
set(gca,'xtick',1:length(sub_n),  'XTickLabel', sub_n)

subplot(2,2,2)
hold on
boxchart(sub_id(keep_idx), n_sub_py(keep_idx), 'BoxFaceColor',c_ord(2,:)); 
scatter(sub_id(keep_idx)+MS_randn_range(1,length(sub_id(keep_idx)),-.1, .1), n_sub_py(keep_idx), 50, 'markerfacecolor', c_ord(2,:), 'markeredgecolor', c_ord(2,:))
set(gca, 'XTickLabel', sub_n)
ylim([0 max([n_Ca1_py(keep_idx) n_sub_py(keep_idx)])])
set(gca,'xtick',1:length(sub_n),  'XTickLabel', sub_n)

subplot(2,2,3)
hold on
boxchart(sub_id(keep_idx), n_Ca1_int(keep_idx), 'BoxFaceColor',c_ord(1,:)); 
scatter(sub_id(keep_idx)+MS_randn_range(1,length(sub_id(keep_idx)),-.1, .1), n_Ca1_int(keep_idx), 50, 'markerfacecolor', c_ord(1,:), 'markeredgecolor', c_ord(1,:))
ylim([0 max([n_Ca1_int(keep_idx) n_sub_int(keep_idx)])])
ylabel('N Fast-spiking Cells')
set(gca,'xtick',1:length(sub_n),  'XTickLabel', sub_n)

subplot(2,2,4)
hold on
boxchart(sub_id(keep_idx), n_sub_int(keep_idx), 'BoxFaceColor',c_ord(2,:)); 
scatter(sub_id(keep_idx)+MS_randn_range(1,length(sub_id(keep_idx)),-.1, .1), n_sub_int(keep_idx), 50, 'markerfacecolor', c_ord(2,:), 'markeredgecolor', c_ord(2,:))
set(gca, 'XTickLabel', sub_n)
ylim([0 max([n_Ca1_int(keep_idx) n_sub_int(keep_idx)])])
set(gca,'xtick',1:length(sub_n),  'XTickLabel', sub_n)