function [psd] = MS_tsd_psd(tsd_in, chan_idx, hann_win, plt_flag)
%% quickly plots a PSD for a channel in a tsd.


if nargin < 2
    chan_idx = size(tsd_in.data,1);
    hann_win = 2^14;
    plt_flag = 1;
elseif nargin < 3
    hann_win = 2^14;
    plt_flag = 1;
elseif nargin < 4
    plt_flag = 1;
end

%% get the psd
for ii = chan_idx
    [psd.(tsd_in.label{ii}(1:end-4)).pxx, psd.(tsd_in.label{ii}(1:end-4)).f] = pwelch(tsd_in.data(ii,:), hanning(hann_win), hann_win/2, hann_win*2 , tsd_in.cfg.hdr{1}.SamplingFrequency);
end

%% plot
c_ord = MS_linspecer(length(chan_idx));

if plt_flag

    figure;
    for ii = chan_idx
        plot(psd.(tsd_in.label{ii}(1:end-4)).f, 10*log10(psd.(tsd_in.label{ii}(1:end-4)).pxx), 'color', c_ord(ii,:));
        hold on;
    end
    xlim([0 300])
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    title('Power Spectral Density');
    legend(tsd_in.label(chan_idx));
    hold off;
end

end