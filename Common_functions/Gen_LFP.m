function data_out = Gen_LFP(cfg_in)
%% Gen_LFP: uses sine waves to generate an artificial LFP signal

% Defaults
% cfg_in.type         = 'sine';   % can be "gaus", 'sine', or 'complex'.  default 'sine';
% cfg_in.nChan        = 4;
% cfg_in.freq         = 60;
% cfg_in.amp          = 1;
% cfg_in.Fs           = 2000;
% cfg_in.len          = 5;
% cfg_in.phase_off    = 90;
% cfg_in.sd           = 6;         % if using the gaussian envolope
% cfg_in.noise        = 'on';
% cfg_in.noise_val    = 0.5;
% cfg_in.complex_freq = [4 34 88 130];
% cfg_in.complex_amp  = [1 1.2 0.8 0.5];
% cfg_in.output       = 'tsd';        % options are tsd, FT, AMPX



%% set defaults
cfg_def = [];
cfg_def.type = 'sine';   % can be "gaus", 'sine', or 'complex'.  default 'sine';
cfg_def.nChan = 4;
cfg_def.freq = 60;
cfg_def.amp = 1;
cfg_def.Fs = 2000;
cfg_def.len = 5;
cfg_def.phase_off = 90;
cfg_def.sd = 6;         % if using the gaussian envolope
cfg_def.noise = 'on';
cfg_def.noise_val = .5;
cfg_def.complex_freq = [4 34 88 130];
cfg_def.complex_amp = [1 1.2 0.8 0.5];
cfg_def.output = 'tsd';        % options are tsd, FT, AMPX
cfg = ProcessConfig2(cfg_def, cfg_in);

%% create the signal
d.tvec = 0:1/cfg.Fs:cfg.len;

switch cfg.type
    
    case 'sine'
        for iChan = 1:cfg.nChan
            if iChan ==1
                d.data(iChan,:) = sin(2*pi*cfg.freq*d.tvec).*cfg.amp;
            else
                if numel(cfg.phase_off) ==1
                    d.data(iChan,:) = sin(2*pi*cfg.freq*d.tvec+degtorad(cfg.phase_off*iChan)).*cfg.amp*iChan;
                else
                    d.data(iChan,:) = sin(2*pi*cfg.freq*d.tvec+degtorad(cfg.phase_off(iChan))).*cfg.amp*iChan;
                end
            end
            if strcmp(cfg.noise, 'on')
                d.data(iChan,:) = d.data(iChan,:)+cfg.noise_val*randn(size(d.tvec));
            end
        end
        
        
    case 'gaus'
        for iChan = 1:cfg.nChan
            if iChan ==1
                d.data(iChan,:) = sin(2*pi*cfg.freq*d.tvec).*cfg.amp*iChan;
                d.data(iChan,:) = d.data(iChan,:).*gausswin(length(d.data(1,:)),cfg.sd)';
            else
                if numel(cfg.phase_off) ==1
                    d.data(iChan,:) = sin(2*pi*cfg.freq*d.tvec+degtorad(cfg.phase_off*iChan)).*cfg.amp*iChan;
                else
                    d.data(iChan,:) = sin(2*pi*cfg.freq*d.tvec+degtorad(cfg.phase_off(iChan))).*cfg.amp*iChan;
                    
                end
                d.data(iChan,:) = d.data(iChan,:).*gausswin(length(d.data(1,:)),cfg.sd)';
            end
            if strcmp(cfg.noise, 'on')
                d.data(iChan,:) = d.data(iChan,:)+cfg.noise_val*randn(size(d.tvec));
            end
        end
        
        
    case 'complex'
        for iChan = 1:cfg.nChan
            freqs =cfg.complex_freq;
            for ifreq =1:length(freqs)
                if iChan ==1
                    temp(ifreq,:) = sin(2*pi*freqs(ifreq)*d.tvec).*cfg.complex_amp(ifreq);
                else
                    if numel(cfg.phase_off) ==1
                        temp(ifreq,:) = (sin(2*pi*freqs(ifreq)*d.tvec+degtorad(cfg.phase_off*iChan)).*cfg.complex_amp(ifreq));
                    else
                        temp(ifreq,:) = (sin(2*pi*freqs(ifreq)*d.tvec+degtorad(cfg.phase_off(iChan))).*cfg.complex_amp(ifreq));
                    end
                end
            end
            d.data(iChan,:) = sum(temp);
            if strcmp(cfg.noise, 'on')
                d.data(iChan,:) = d.data(iChan,:)+cfg.noise_val*randn(size(d.tvec));
            end
        end
end
% make labels
for iChan = 1:cfg.nChan
    d.label{iChan} = ['CSC' num2str(iChan)];
end
%% convert to output format

switch cfg.output
    
    case 'tsd' % time series data
        d.cfg = cfg;
        for nChan = cfg.nChan:-1:1
            d.cfg.hdr{nChan}.SamplingFrequency = cfg.Fs;
            d.cfg.hdr{nChan}.FileType = 'CSC';
            d.cfg.hdr{nChan}.AcqEntName = ['CSC' num2str(nChan)];
        end
        d.type = cfg.output;
        d.cfg.SessionID = 'ARTIFICIAL';
        d.cfg.history.mfun = []; d.cfg.history.cfg = [];
        data_out = d;
        figure(9999)
        subplot(2,1,1)
        plot(data_out.tvec, data_out.data)
        subplot(2,1,2)
        plot(data_out.tvec, data_out.data)
        xlim([0 cfg.len/25])
        legend(data_out.label, 'location', 'south', 'orientation', 'horizontal')
        
    case 'FT' % fieldtrip
        data_out = TSDtoFT([],d);
        figure(9999)
        subplot(2,1,1)
        plot(data_out.time{1}, data_out.trial{1})
        subplot(2,1,2)
        plot(data_out.time{1}, data_out.trial{1})
        xlim([0 cfg.len/25])
end
%% plot check

