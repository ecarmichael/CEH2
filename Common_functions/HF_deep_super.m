function [S_out, rip_out] = HF_deep_super(S, csc, swr_iv, probe, chan_idx)

%% HF_deep_super: uses the SWR times to create a ripple average for all the
%  channels on a shank by shank basis. Determines deep vs superficial 
%  cutoff using the methods from Petersen et al. 2022 with the addition of
%  a Hampel outlier identifier (1 sample) and fitting a sine wave to the
%  mean voltage in the -40-12ms range.  
%
%
%
%    Inputs:
%    - S [struct]  spike times in the TS format
%
%    - csc [struct]  LFP in the tsd format
%
%    - swr_iv [struct]  swr events in the iv format
%
%    - probe [string]  can be either 'A4x16', A5x12', 'Buz32'
%
%    - chan_idx [n x 1] indices of channels to use.
%
%    Outputs:
%    - S_out [struct]  spike TS with the deep/super classification
%
%
%
%
% EC 2026-09-01   initial version
%
%
%
%% initialize

if nargin < 5
    chan_idx = size(csc.data,1);
end


%% setup the probe
if strcmpi(probe, 'A4x16')
    shank{1} = 1:16;
    shank{2} = 17:32;
    shank{3} = 33:48;
    shank{4} = 49:64;
    ycoords   = repmat([100, 50, 0:-20:-260],1,4);
    ycoords   = ycoords(:);
    x_space = [8.66 8.66, repmat([0 17.3], 1,7)]; % space between probes x in um;
    xcoords   = [x_space, x_space , x_space , x_space];                    %repmat([1 2 3 4]', 1, Nchannels/4);
    xcoords   = xcoords(:);
    ref_idx = [1 2 17 18 33 34 49 50]; 
elseif strcmpi(probe, 'Buz32')
    shank{1} = 1:8;
    shank{2} = 9:16;
    shank{3} = 17:24;
    shank{4} = 25:32;
    ycoords   = repmat([0:-40:-120, -140 -100:40:-20],1,4);
    ycoords   = ycoords(:);
    xcoords = repmat([0 8.5 17 17+8.5 17+8.5*2 34+8.5 34+17 51+8.5 ],1,4); % space between probes x in um;
    xcoords   = xcoords(:);
ref_idx = []; 
elseif strcmpi(probe, 'A5x12')
    shank{1} = 1:12;
    shank{2} = 13:24;
    shank{3} = 24:40;
    shank{4} = 41:52;
    shank{5} = 53:64;
    ycoords = [repmat(0:-10:-110,1,2), [1500 1000 500 0:-100:-1200], repmat(0:-10:-110,1,2)];
    ycoords = ycoords(:);
    x_space = repmat([0 20], 1,6); % space between probes x in um;
    xcoords   = [x_space, x_space , (ones(1,16))+10, x_space , x_space ];                    %repmat([1 2 3 4]', 1, Nchannels/4);
    xcoords   = xcoords(:);
    ref_idx = []; 
end

n = length(shank);

% get a subset of ripples
[~, idx] = sort(swr_iv.usr.mean_filt, 'descend');
if length(swr_iv.tstart) < 100

    idx = 1:length(swr_iv.tstart);
else
    idx = idx(10:100);
end
r_iv = SelectIV([], swr_iv, idx);

rip_cent = IVcenters(r_iv);

%% filter the LFP into the ripple band.

cfg_swr.type = 'cheby1'; %Cheby1 is sharper than butter
cfg_swr.f  = [125 200]; % broad, could use 150-200?
cfg_swr.order = 4; %type filter order (fine for this f range)
cfg_swr.display_filter = 0; % use this to see the fvtool

csc_f = FilterLFP(cfg_swr, csc);
%% loop over shanks

% subplot index
% s_idx = [1 2 3; 5 6 7; 9 10 11; 13 14 15; 17 18 19];
% p_idx  = 4:4:20;
c_ord = MS_linspecer(length(shank{1}));

% flip_idx = []; 
for ii = 1:length(shank)

    % channels for the current shank
    this_ch = chan_idx(shank{ii});
    rm_idx = find(ismember(this_ch, ref_idx)); %Remove ref channels if needed.
    
    c_ord(rm_idx,:) = repmat([.5 .5 .5], length(rm_idx), 1);

    win = .125;
    rippleAvg = [];
    ripplePow = [];
    rippleAvg_f = [];
    % get the average
    for jj = length(rip_cent):-1:1

        % get the ripple power;
        rip_idx = nearest_idx([r_iv.tstart(jj), r_iv.tend(jj)],csc.tvec);

        ripplePow(:,jj) = mean(abs(hilbert(csc_f.data(this_ch, rip_idx(1):rip_idx(2)))).^2,2);


        win_idx = nearest_idx([rip_cent(jj)-win, rip_cent(jj)+win],csc.tvec);
        rippleAvg(:, :, jj) = csc.data(this_ch, win_idx(1):win_idx(2));
        rippleAvg_f(:, :, jj) = csc_f.data(this_ch, win_idx(1):win_idx(2));

    end

    % collect and average
    tvec = csc.tvec(win_idx(1):win_idx(2));
    tvec = tvec - tvec(1);
    tvec = tvec - win;

    rippleAvg = mean(rippleAvg,3);
    rippleAvg_f = mean(rippleAvg_f,3);

    
    % get the flip point based on Petersen et al. The polarity of the 
    % sharp-wave is determined from an interval before the average ripple 
    % peak at -40.8ms to -12.8ms, as this has shown to be rebust across 
    % species (rats and mice). The algorith will look for the point where 
    % the polarity of the sharp-wave flips.

    d_idx = nearest_idx([-.030, -0], tvec);
    ripple_diff = mean(rippleAvg(:,d_idx(1):d_idx(2)),2); 
    this_idx = ismember(1:length(ripple_diff), rm_idx); 

    % apply hampel outlier detector for bad channels
    y = hampel(ripple_diff(~this_idx),1);


    % fit a sine wave for smoothing
    fit_out = fit((1:length(ripple_diff(~this_idx)))', y, 'sin1');
    s_fit= fit_out((1:length(ripple_diff(~this_idx))));



    % find the crossing point (first above 0)

    flip_idx = find(diff(sign(s_fit)) ~= 0);
    if isempty(flip_idx) && sum(sign(s_fit) >0 ) == length(s_fit) % if everything above the later
        deep_idx = true(size(this_ch));
        super_idx = logical(~deep_idx);
    elseif isempty(flip_idx) && sum(sign(s_fit) < 0 ) == length(s_fit) % if everything is below. 
        deep_idx = false(size(this_ch));
        super_idx = logical(~deep_idx);
    else
        deep_idx = zeros(size(this_ch));
        deep_idx(1:flip_idx+length(rm_idx)) = 1;

        deep_idx = logical(deep_idx);
        super_idx = logical(~deep_idx);
    end
    s_fit = [NaN NaN s_fit'];


    flip_ch = find(deep_idx); 
    flip_ch = flip_ch(end); 

    figure(10+ii)
    clf

    % raw Ripple triggered average with spacing. 
    subplot(2,3,[1 4])
    cla;
    hold on
    for kk = 1:length(this_ch)
        if kk == flip_ch
            plot(tvec, (rippleAvg(kk, :)*.1) + ycoords(this_ch(kk)), 'color', 'k', 'linewidth', 2)
        else
            plot(tvec, (rippleAvg(kk, :)*.1) + ycoords(this_ch(kk)), 'color', c_ord(kk,:))
        end
    end
    ylim([min(ycoords(this_ch))-50 max(ycoords(this_ch))+50])
    set(gca, "XTick", -.12:.04:.12, 'XtickLabel', [-.12:.04:.12]*1000)


    % no offset a la Mizuseki 2014
    subplot(2,3,2)
    cla;
    hold on
    for kk = 1:length(this_ch)
        if kk == flip_ch
            plot(tvec, rippleAvg(kk, :), 'color', 'k')
        else
            plot(tvec, rippleAvg(kk, :), 'color', c_ord(kk,:))
        end
    end
    xlim([-.080 .08])
    set(gca, "XTick", -.08:.04:.08, 'XtickLabel', [-.08:.04:.08]*1000)
    xline([-.030, -0])

    % make the probe out of rectangles
    fac = 800;
    for kk = 1:length(this_ch)

        x_off = (xcoords(this_ch(kk))/fac) -.07-min(xcoords(this_ch)/fac);
        rectangle('Position', [x_off , (ycoords(this_ch(kk)))-300, 6/fac, 12], ...
            'FaceColor', c_ord(kk,:), 'EdgeColor', 'none');

        if kk == flip_ch
            text(x_off, (ycoords(this_ch(kk)))-300, num2str(this_ch(kk)), VerticalAlignment='middle', HorizontalAlignment='right', fontweight = 'bold')
        else
            text(x_off, (ycoords(this_ch(kk)))-300, num2str(this_ch(kk)), VerticalAlignment='middle', HorizontalAlignment='right')
        end
    end

    % add the filtered means
    subplot(2,3,5)
    cla;
    hold on
    for kk = 1:length(this_ch)
        plot(tvec, (rippleAvg_f(kk, :)*.5) + ycoords(this_ch(kk)), 'color', c_ord(kk,:))
    end
    ylim([min(ycoords(this_ch))-50 max(ycoords(this_ch))+50])



    subplot(2,3, 3)
    cla
    hold on
    b=bar(ycoords(this_ch), ripple_diff');
    b.FaceColor = 'flat';
    b.CData = flipud(c_ord);
    b.EdgeColor = 'none';
    
    plot(ycoords(this_ch), s_fit, '-k')
    scatter(ycoords(deep_idx), s_fit(deep_idx), 50, 'filled','markerfacecolor',  'b')
    scatter(ycoords(super_idx), s_fit(super_idx), 50, 'filled','markerfacecolor',  'r')


    % eb = errorbar(ycoords(this_ch), mean(ripplePow,2, 'omitnan'), MS_SEM_vec(ripplePow'), 'vertical');
    % eb.LineStyle = 'none';
    % eb.Color = [.2 .2 .2];
    % eb.LineWidth =1;


    view(90,90)
    xlim([min(ycoords(this_ch))-50 max(ycoords(this_ch))+50])
    set(gca, 'xDir', 'reverse')

    

    subplot(2,3, 6)
    cla
    hold on
    b=bar(ycoords(this_ch), mean(ripplePow, 2, 'omitnan')');
    b.FaceColor = 'flat';
    b.CData = flipud(c_ord);
    b.EdgeColor = 'none';

    eb = errorbar(ycoords(this_ch), mean(ripplePow,2, 'omitnan'), MS_SEM_vec(ripplePow'), 'vertical');
    eb.LineStyle = 'none';
    eb.Color = [.2 .2 .2];
    eb.LineWidth =1;


    view(90,90)
    xlim([min(ycoords(this_ch))-50 max(ycoords(this_ch))+50])
    ylim([0 inf])
    set(gca, 'xDir', 'reverse')

% collect the outputs
rip_out{ii} = []; 
rip_out{ii}.deep = deep_idx; 
rip_out{ii}.super = super_idx; 
rip_out{ii}.flip = flip_idx;

rip_out{ii}.rippleAvg = rippleAvg;
rip_out{ii}.rippleAvg_f = rippleAvg_f;
rip_out{ii}.ripplePow = ripplePow; 
rip_out{ii}.ycoords = ycoords(this_ch); 

end


%% apply the cutoffs to to the spike file using the channel with the maximum spike amplitude. 
deep_chan = []; 

for ii = 1:length(rip_out)
    this_ch = chan_idx(shank{ii});

    deep_chan = [deep_chan this_ch(rip_out{ii}.deep)]; 

end

S_out = S; 

S_out.usr.deep = ismember(S_out.usr.ch, deep_chan); 