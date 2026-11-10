function CSDs = HF_probe_CSD(data_in, band, evts)
%% HF_probe_CSD:  compute a shank by shank csd using one linear set of electrodes per shank.
%
%
%
%    Inputs:
%    - data_in [struct]  needs to contain the CSC and a set of probe names.
%
%
%
%
%    Outputs:
%    -
%
%
%
%
% EC 2026-11-09   initial version
%
%
%
%% initialize

if nargin < 2
    band = 'SWR';
    evts = [];
elseif nargin < 3
    evts = [];
end

%% define the probes using the common pairings

% if strcmpi(data_in.probe{1}, '4x16')  & strcmpi(data_in.probes, 'Buz32')
A4x16.shank{1} = 1:16;
A4x16.shank{2} = 17:32;
A4x16.shank{3} = 33:48;
A4x16.shank{4} = 49:64;
ycoords   = repmat([100, 50, 0:-20:-260],1,4);
A4x16.ycoords   = ycoords(:);
x_space = [8.66 8.66, repmat([0 17.3], 1,7)]; % space between probes x in um;
xcoords   = [x_space, x_space , x_space , x_space];                    %repmat([1 2 3 4]', 1, Nchannels/4);
A4x16.xcoords   = xcoords(:);
xcoords_off = [x_space, x_space+150 , x_space+300 , x_space+450]; % with the offsets.
A4x16.xcoords_off = xcoords_off(:);
A4x16.ref_idx = [1 2 17 18 33 34 49 50];
% buz 32

buz32.shank{1} = 1:8;
buz32.shank{2} = 9:16;
buz32.shank{3} = 17:24;
buz32.shank{4} = 25:32;
ycoords   = repmat([0:-40:-120, -140 -100:40:-20],1,4);
buz32.ycoords   = ycoords(:);
xcoords = repmat([0 8.5 17 17+8.5 17+8.5*2 34+8.5 34+17 51+8.5 ],1,4); % space between probes x in um;
buz32.xcoords   = xcoords(:);
xcoords_off = [[0 8.5 17 17+8.5 17+8.5*2 34+8.5 34+17 51+8.5 ], [0 8.5 17 17+8.5 17+8.5*2 34+8.5 34+17 51+8.5 ]+200,...
    [0 8.5 17 17+8.5 17+8.5*2 34+8.5 34+17 51+8.5 ]+400 , [0 8.5 17 17+8.5 17+8.5*2 34+8.5 34+17 51+8.5 ]+600]; % with the offsets.
buz32.xcoords_off = xcoords_off(:);
buz32.ref_idx = [];
% elseif strcmpi(data_in.probes, '5x12')
A5x12.shank{1} = 1:16;
A5x12.shank{2} = 17:28;
A5x12.shank{3} = 29:40;
A5x12.shank{4} = 41:52;
A5x12.shank{5} = 53:64;
ycoords = [ [1500 1000 500 0:-100:-1200], repmat(0:-10:-110,1,2), repmat(0:-10:-110,1,2)];
A5x12.ycoords = ycoords(:);
x_space = repmat([0 20], 1,6); % space between probes x in um;
xcoords   = [(ones(1,16))+10, x_space, x_space, x_space , x_space ];                    %repmat([1 2 3 4]', 1, Nchannels/4);
A5x12.xcoords   = xcoords(:);
xcoords_off   = [ (ones(1,16)*320)+10, x_space, x_space+160, x_space+480, x_space+640];                    %repmat([1 2 3 4]', 1, Nchannels/4);
A5x12.xcoords_off = xcoords_off(:);
A5x12.ref_idx = [1 2 3];


% set the probes based on config
if strcmpi(data_in.probe{1}, '5x12')  &  strcmpi(data_in.probe{2}, 'buz32')

    probe.shank = A5x12.shank;
    for ii = 1:length(buz32.shank)
        probe.shank{end+1} = buz32.shank{ii}+64;
    end
    probe.xcoords = [A5x12.xcoords; buz32.xcoords];
    probe.ycoords = [A5x12.ycoords; buz32.ycoords];
    probe.xcoords_off = [A5x12.xcoords_off; buz32.xcoords_off];
    probe.ref_idx = [A5x12.ref_idx; buz32.ref_idx+64];

elseif strcmpi(data_in.probe{1}, '5x12')  &  strcmpi(data_in.probe{2}, '4x16')

    probe.shank = A5x12.shank;
    for ii = 1:length(A4x16.shank)
        probe.shank{end+1} = A4x16.shank{ii}+64;
    end
    probe.xcoords = [A5x12.xcoords; A4x16.xcoords];
    probe.ycoords = [A5x12.ycoords; A4x16.ycoords];
    probe.xcoords_off = [A5x12.xcoords_off; A4x16.xcoords_off];
    probe.ref_idx = [A5x12.ref_idx; A4x16.ref_idx+64];

elseif strcmpi(data_in.probe{1}, '4x16')  &  strcmpi(data_in.probe{2}, 'buz32')

    probe.shank = A5x12.shank;
    for ii = 1:length(A4x16.shank)
        probe.shank{end+1} = A4x16.shank{ii}+64;
    end
    probe.xcoords = [A5x12.xcoords; A4x16.xcoords];
    probe.ycoords = [A5x12.ycoords; A4x16.ycoords];
    probe.xcoords_off = [A5x12.xcoords_off; A4x16.xcoords_off];
    probe.ref_idx = [A5x12.ref_idx; A4x16.ref_idx+64];

else
    error('Unknown probe configuration')
end

clearvars -except data_in probe evts band

%% loop over the probe shanks and get the CSD


% for iS = 1:length(probe.shank)
iS = 1

this_idx = probe.shank{iS};

%remove reference electrodes
this_idx(probe.ref_idx) = [];

% grab the channels of interest
this_csc = MS_select_tsd(data_in.csc, this_idx);

% get the positions (y only).
this_pos = probe.ycoords(this_idx);

if contains(lower(band), 'swr')
    % filer into band of interest
    cfg_filter.type = 'cheby1'; %Cheby1 is sharper than butter
    cfg_filter.f  = [125 200]; % broad, could use 150-200?
    cfg_filter.order = 4; %type filter order (fine for this f range)
    cfg_filter.display_filter = 0; % use this to see the fvtool

    csc_f = FilterLFP(cfg_filter, this_csc);

    d_fs = floor(1/(50/csc_f.cfg.hdr{1}.SamplingFrequency)); % 50ms window. 
    data = smoothdata(abs(hilbert(csc_f.data)),2, "gaussian", d_fs); 

elseif contains(lower(band), 'low gamma')

    cfg_filter = [];
    cfg_filter.f = [30 55];
    cfg_filter.type = 'butter';
    cfg_filter.order = 3;

    csc_f = FilterLFP(cfg_filter, this_csc);
    data = csc_f.data; 

elseif contains(lower(band), 'mid gamma')

    cfg_filter = [];
    cfg_filter.f = [65 90];
    cfg_filter.type = 'butter';
    cfg_filter.order = 3;

    csc_f = FilterLFP(cfg_filter, this_csc);
    
    data =csc_f.data;

elseif contains(lower(band), 'high gamma')

    cfg_filter = [];
    cfg_filter.f = [100 250];
    cfg_filter.type = 'butter';
    cfg_filter.order = 3;

    csc_f = FilterLFP(cfg_filter, this_csc);

    data =csc_f.data;

elseif contains(lower(band), 'theta')

    cfg_filter = [];
    cfg_filter.f = [6 12];
    cfg_filter.type = 'cheby1';
    cfg_filter.order = 3;
    csc_f = FilterLFP(cfg_filter, this_csc);
    data = csc_f.data; 
end
    % compute CSD
    [csd, csd_inds]= csd5pt_Schomburg(data, abs(mode(diff(this_pos))));


    % smooth the 
    % tsd the CSD for easier restriction
    csd_tsd{iS} = tsd(csc_f.tvec, csd,num2cell(csd_inds));
    % end


    %% get the mean across events

    if ~isempty(evts.tstart)
        win = .5;
        win_idx = win*this_csc.cfg.hdr{1}.SamplingFrequency;
        evt_cent = IVcenters(evts);
        evt_idx = find(ismember(csd_tsd{iS}.tvec, evt_cent));

        evt_tvec = -win:1/this_csc.cfg.hdr{1}.SamplingFrequency:win;

        %% average the csd across events
        all_csd = []; all_evt_csc =[];  all_csd_peak = []; 
        for ee = length(evt_idx):-1:1

            all_csd(:,:,ee) = csd_tsd{iS}.data(:,evt_idx(ee)-win_idx:evt_idx(ee)+win_idx);
            all_csd_peak(:,ee) = csd_tsd{iS}.data(:,evt_idx(ee));

            % lds_method
            all_evt_csc(:,:,ee) = data(:,evt_idx(ee)-win_idx:evt_idx(ee)+win_idx);

        end

        % LDS 2026 method
        mean_evt_csc = mean(all_evt_csc,3);
        lds_csd = NaN(size(mean_evt_csc,1)-2, size(mean_evt_csc,2));
        for iC = 2:size(mean_evt_csc,1)-1
            for jj = 1:length(mean_evt_csc)
                lds_csd(iC, jj) = -(mean_evt_csc(iC-1,jj)-2*mean_evt_csc(iC,jj)+mean_evt_csc(iC+1,jj));
            end
        end
        lds_csd(1,:) = []; 

        % interp across channels
        mean_evt_csc_int =[]; 
        int_pos = fliplr(this_pos(end):25:this_pos(1))'; 
        % og_pos = this_pos(2:size(mean_evt_csc,1)-1); 
        % int_pos =  fliplr(og_pos(end):20:og_pos(1));
        for jj = length(mean_evt_csc):-1:1
            mean_evt_csc_int(:,jj) = interp1(mean_evt_csc(:,jj), this_pos, int_pos);
        end

        % lds_csd_int = []; 
        % 
        % for jj = length(mean_evt_csc):-1:1
        %     lds_csd_int(:,jj) = interp1(lds_csd(:, jj), og_pos, int_pos); 
        % end
        % 
        % lds_csd = smoothdata(lds_csd, 1,"gaussian",50/abs(mode(diff(this_pos)))); % smooth vertically at 50um
        % lds_csd_int = smoothdata(lds_csd_int, 1,"gaussian",50/abs(mode(diff(int_pos)))); % smooth vertically at 50um

        lds_csd_peak(:,ee) = lds_csd(:,floor(win_idx/2));


    end
    %% check plot

    csc_r = restrict(this_csc, 0 ,1500);
    csc_f_r = restrict(csc_f, 0 ,1500);
    csd_r = restrict(csd_tsd{iS}, 0, 1500);

    figure(100+iS)
    clf;
    ax(1) = subplot(2,4,1:2);

    imagesc(csd_r.tvec,csd_inds, csd_r.data)
    set(gca, 'YTickLabel', this_pos(csd_inds), 'YDir', 'reverse')

    ax(2) = subplot(2,4,3:4);
    hold on
    lab =[]; y_t = [];
    for ii = 1:size(csc_r.data,1)


        plot(csc_r.tvec, csc_r.data(ii,:)/1000+ii);
        plot(csc_f_r.tvec, csc_f_r.data(ii,:)/1000+ii, 'k');

        lab{ii} = csc_r.label{ii};
        y_t(ii) = ii;
    end
    set(gca, 'YTick', y_t, 'YTickLabel', this_pos, 'YDir', 'reverse')

    linkaxes(ax, 'xy')

    % ylim([1:size(csc_r.data,1)])

    if exist('all_csd', 'var')
        subplot(2,4,5)
        cla
        hold on
        imagesc(evt_tvec, csd_inds,mean(all_csd, 3));
        % add in the signals. 
        for ii = 1:size(mean_evt_csc,1)
            plot(evt_tvec, (mean_evt_csc(ii,:)/50)+ii);
        end
        set(gca, 'YTick', y_t, 'YTickLabel', this_pos, 'YDir', 'reverse')

        xlabel('Time (s)');
        ylabel('CSD position');
        title('Event-averaged CSD');
        xlim([-.1 .1])

        subplot(2,4,6)
        cla
        plot( normalize(mean(all_csd_peak,2)), csd_inds, '-.');
        set(gca, 'YTickLabel', this_pos(csd_inds), 'YDir', 'reverse')
        ylabel('CSD position');
        title('CSD at peak');
        xlim([-1 1])

        % same thing with the LDS version
        subplot(2,4,7)
        cla
        hold on
        imagesc(evt_tvec, 2:length(this_pos)-1,lds_csd);
        % add in the signals. 
        for ii = 1:size(mean_evt_csc,1)
            plot(evt_tvec, (mean_evt_csc(ii,:)/50)+ii);
        end
        set(gca, 'YTick', y_t, 'YTickLabel', this_pos, 'YDir', 'reverse')
        xlabel('Time (s)');
        ylabel('CSD position');
        title('LDS CSD');
        xlim([-.1 .1])

        subplot(2,4,8)
        cla
        plot( normalize(lds_csd_peak), 2:length(this_pos)-1, '-.');
        set(gca, 'YTickLabel', this_pos(2:length(this_pos)-1), 'YDir', 'reverse')
        ylabel('CSD position');
        title('LDS CSD peak');
        xlim([-1 1])

    end
end
