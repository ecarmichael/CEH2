function [A_Temp, time_proj] = MS_PCA_ICA(cfg_in, S, pos, csc)
%% MS_PCA_ICA:
%
%
%
%    Inputs:
%    -
%
%
%
%    Outputs:
%    -
%
%
%
%
% EC 2023-09-10   initial version
%
%
%
%% initialize

if nargin <3
    pos = [];
    csc= [];
elseif nargin < 4
    csc = [];
end


cfg_def = [];
cfg_def.pos_ass = 1; % restrict only to assemblies with significantly positive weights
cfg_def.pos_ass_z = 1.98; % z score threshold for positive weights.
cfg_def.bin_s = .001;
cfg_def.mov = 0;
cfg_def.spd_th = 2.5; % movement threshold in cm/s;
cfg_def.plot = 0;
cfg_def.lfp_idx = 1; % just plot the first lfp channel.
cfg_def.t_win = []; % fill this in if you want a time course without a using pos.

cfg = ProcessConfig(cfg_def, cfg_in);


%%

if ~isempty(pos)

    S = restrict(S, pos.tvec(1), pos.tvec(end));
    if ~isempty(csc)
        csc = restrict(csc, pos.tvec(1), pos.tvec(end));

    end

    speed = pos;

    speed.data = pos.data(end-1,:);


    if cfg.mov
        %     move_idx = speed> cfg.spd_th;
        cfg_mov = []; cfg_mov.method = 'raw'; cfg_mov.operation = '>'; cfg_mov.threshold = cfg.spd_th; % speed limit in cm/sec
        iv_fast = TSDtoIV(cfg_mov,speed); % only keep intervals with speed above thresh

        S = restrict(S, iv_fast);
    end
end


%%
if ~isempty(cfg.t_win)
    tstart = cfg.t_win(1);
    tend = cfg.t_win(end);
else
    tstart = pos.tvec(1);
    tend = pos.tvec(end);
end
tbin_edges = tstart:cfg.bin_s:tend; % vector of time bin edges (for histogram)
tbin_centers = tbin_edges(1:end-1)+cfg.bin_s/2; % vector of time bin centers (for plotting)

data_h = [];
for ii = size(S.t,2):-1:1


    spk_count = histc(S.t{ii},tbin_edges); % get spike counts for each bin
    spk_count = spk_count(1:end-1); % ignore spikes falling exactly on edge of last bin.
    spk_count = smoothdata(spk_count, 10, 'gaussian');

    data_h(:,ii) = spk_count;
end

% figure
% hold on
% for ii = 1:size(data_h,2)
% plot(data_h(:,ii)+ii*10)
%
% end
%%
rng default
A_Temp = assembly_patterns(data_h');

if isempty(A_Temp)
    A_temp= [];
    time_proj = [];
    return
end

rng default
time_proj = assembly_activity(A_Temp,data_h');



%% get the positive weight cells

for ii = 1:size(A_Temp,2)

    if cfg.pos_ass
        A_pos_idx{ii} = find(zscore(A_Temp(:,ii)) > 1.98);
    else
        A_pos_idx{ii} = 1:length(A_Temp(:,ii));
    end
end


if cfg.pos_ass
    keep_idx = cellfun(@length, A_pos_idx) > 0;
    fprintf('<strong>%s</strong>:<strong> %0.f / %0.f</strong> Assemblies detected had significantly positive weights.\n', mfilename, sum(keep_idx), length(keep_idx))

else
    keep_idx = ones(size(A_pos_idx));
    fprintf('<strong>%s</strong>:<strong> %0.f</strong> Assemblies detected \n', mfilename,  length(keep_idx))

end



A_idx = find(keep_idx);
A_pos_idx(~keep_idx) = [];
A_pos = A_Temp;
A_pos(:,~keep_idx) = [];
time_proj_pos = time_proj;
time_proj_pos(~keep_idx,:) = [];

%% plot the positive assemblies;

if cfg.plot
    figure(104);clf; hold on

    c_ord = parula(size(A_pos,2)+2);


    for ii = size(A_pos,2):-1:1
        figure(104)
        subplot( ceil(size(A_pos,2)/4),4,ii)
        hold on
        stem(A_pos(:,ii), 'color', c_ord(ii,:))
        view(90,90)

        stem(A_pos_idx{ii}, A_pos(A_pos_idx{ii},ii), 'color', c_ord(ii,:), 'MarkerFaceColor', c_ord(ii,:))


        title(['Assembly#' num2str(A_idx(ii))])
    end

end

%% sample plot
% time_proj_tsd = tsd(tbin_centers, time_proj_pos);
%
% cfg_plot = [];
% cfg_plot.spkColor = parula(length(S.t));
% cfg_plot.LineWidth = 2;
% cfg_plot.lfp = time_proj_tsd;

%%

if cfg.plot
    figure(900)
    clf
    maximize
    if ~isempty(pos)

        ax(1) = subplot(7, 1, 1);

    end
    hold on
    plot(pos.tvec,pos.data(1:2,:));
    xlim([pos.tvec(1) pos.tvec(end)])
    ylabel('position on track (cm)')
    set(gca, 'XTick', []);

    if ~isempty(pos)
        set(gca, 'xtick', [])
        ax(2) = subplot(7,1,2:5);
    else
        ax(2) = subplot(7,1,1:5);
    end
    cla
    hold on
    off_set = 0; these_idx = [];
    for iA = 1:length(A_pos_idx)

        for ii = size(A_pos_idx{iA}, 1):-1:1
            iC = A_pos_idx{iA}(ii);
            plot([S.t{iC}, S.t{iC}]', [(ones(size(S.t{iC}))*ii)-.5+off_set, (ones(size(S.t{iC}))*ii)+.5+off_set]', 'color', c_ord(iA,:), 'linewidth', 2)
        end
        off_set = off_set+size(A_pos_idx{iA}, 1);
        these_idx = [these_idx, A_pos_idx{iA}']; % keep track to avoid overlap;
    end
    non_ass_idx = 1:length(S.t);
    rm_idx = (ismember(non_ass_idx, these_idx));
    non_ass_idx(rm_idx) = [];

    for ii = size(non_ass_idx, 2):-1:1
        if ~isempty(S.t{non_ass_idx(ii)})
            plot([S.t{non_ass_idx(ii)}, S.t{non_ass_idx(ii)}]', [(ones(size(S.t{non_ass_idx(ii)}))*ii)-.5+off_set, (ones(size(S.t{non_ass_idx(ii)}))*ii)+.5+off_set]', 'color', [.7 .7 .7 .7], 'linewidth', 2)
        end
    end
    ylim([0 size(S.t, 2)])
    set(gca, 'XTick', [], 'YDir', 'normal', 'ytick', [])
    ylabel('Cell activity')

    if isempty(csc)
        ax(3) = subplot(7,1,6:7);
    else
        ax(3) = subplot(7,1,6);
    end
    cla
    hold on
    leg_val = [];
    for iA = 1:length(A_pos_idx)
        plot(tbin_centers, zscore(time_proj(iA,:)), 'color', c_ord(iA,:))

        leg_val{iA} = ['Assembly #' num2str(iA)];
    end

    ylabel({'zscore'; 'react strength'})
    xlabel('time (s)')

    legend(leg_val, 'Orientation', 'horizontal')


    if ~isempty(csc)
        set(gca, 'xtick', [])
        ax(4) = subplot(7,1,7);
    end
    plot(csc.tvec, csc.data(cfg.lfp_idx,:));


    linkaxes(ax, 'x')
    xlim([pos.tvec(1) pos.tvec(end)])
end


%% plot the assembly spatial maps


for iA = 1:length(A_pos_idx)

    this_S = S;
    this_S.t = [];
    this_S.label = [];
    this_S.t = S.t(A_pos_idx{iA});
    this_S.label = S.label(A_pos_idx{iA});
    cfg_place = [];
    cfg_place.plot = 0;
    TCs{iA} = MS_get_place_field(cfg_place, this_S, pos, speed);



end

if cfg.plot
    figure(99)
    cla
    n = length(A_pos_idx);
    m = max(cellfun(@length, A_pos_idx))+2;
    splt = reshape(1:n*m, m, n)';


    for iA = 1:length(A_pos_idx)
        these_tcs = [];
        for ii = 1:length(A_pos_idx{iA})
            subplot(n, m, splt(iA, ii))

            imagesc(TCs{iA}.tc{ii})
            these_tcs = cat(3,these_tcs, TCs{iA}.tc{ii});
            set(gca, 'XTick', [], 'YTick', [])
            if ii == 1
                ylabel(['Assembly #' num2str(iA)])
            end
        end
        subplot(n, m, splt(iA, m))
        imagesc(mean(these_tcs,3))
        set(gca, 'XTick', [], 'YTick', [])
        title('Assembly mean')
    end



end


end


