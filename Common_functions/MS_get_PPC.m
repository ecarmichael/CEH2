function MS_get_PPC(cfg_in, spike_fname, lfp_fname)
%% MS_get_PPC: wrapper function to compute Pair-wise Phase Consistency using the FieldTrip toolbox.  Will also run shuffle comparison.
%
%
%
%    Inputs:
%    - cfg_in: [struct] configurations
%
%    - spike_fname name of the spike file to process
%
%    - csc [struct] contains LFP data
%
%
%    Outputs:
%    -
%
%
%
%
% EC 2021-07-29   initial version
%
%
%
%% initialize
cfg_def = [];
cfg_def.plot = 0;
cfg_def.shuffle = 100;
cfg_def.inter_dir = [cd filesep 'PPC'];

cfg = ProcessConfig(cfg_def, cfg_in);


if ~exist(cfg.inter_dir)
    mkdir(cfg.inter_dir)
end
%%

spike = ft_read_spike(spike_fname); % needs fixed read_mclust_t.m
if isempty(spike.unit)
    spike.unit{1} = ones(1,length(spike.timestamp{1}));
    if exist([spike_fname(1:6) '-wv.mat'], 'file')
        load([spike_fname(1:6) '-wv.mat'], 'mWV');
        spike.waveform{1} = repmat(mWV,1,1, length(spike.timestamp{1}));
    end
end
fc = lfp_fname;
data = ft_read_neuralynx_interp(fc);
data_all = ft_appendspike([],data, spike);

if sum(isnan(data_all.trial{1}(1,:))) >0
    nan_blocks = diff(isnan(data_all.trial{1}(1,:)));
    
    cfg_trl = [];
    cfg_trl.trl = [data_all.hdr.Fs find(nan_blocks < 0)+1; find(nan_blocks > 0), length(data_all.time{1})]';
    cfg_trl.trl(:,3) = zeros(size(cfg_trl.trl,1),1);
    
    data_trl = ft_redefinetrial(cfg_trl, data_all);
    Blocks = {'Pre_sleep', 'W_maze', 'OF', 'Post_sleep'};
else
    
    cfg_trl = [];
    cfg_trl.begsample = 1;
    cfg_trl.endsample = length(data_all.time{1});
    data_trl = ft_redefinetrial(cfg_trl, data_all);
    Blocks = {'whole session'};
end



%%

spk_chan = data_trl.label{2};
lfp_chan = fc{1}(1:end-4);

cfg_i              = [];
cfg_i.timwin       = [-0.002 0.006]; % remove 4 ms around every spike
cfg_i.spikechannel = spk_chan;
cfg_i.channel      = lfp_chan;
cfg_i.method       = 'linear'; % remove the replaced segment with interpolation

data_i        = ft_spiketriggeredinterpolation(cfg_i, data_trl);

%%
for iB = 1:length(Blocks)
    close all
    cfg_this_trl = [];
    cfg_this_trl.trials = iB;
    
    data_this_trl = ft_redefinetrial(cfg_this_trl, data_i);
    
    
    if isfield(cfg, 'plot')
    figure(iB);
    subplot(221)
    plot(data_this_trl.time{1}(1:200000), data_this_trl.trial{1}(1,1:200000))
    hold on
    plot(data_this_trl.time{1}(1:200000), data_this_trl.trial{1}(2,1:200000)*500)
    xlabel('Sample data from 0-100s');
    end

    %% STA
    cfg_sta              = [];
    cfg_sta.timwin       = [-0.5 0.5]; %
    cfg_sta.spikechannel = spk_chan;
    cfg_sta.channel      = lfp_chan;
    %     cfg_sta.keeptrials   = 'yes';
    PPC.(Blocks{iB}).staAll  = ft_spiketriggeredaverage(cfg_sta, data_this_trl);
    
    % plot
    if isfield(cfg, 'plot')
        
        figure(iB)        
        subplot(223)
        plot(PPC.(Blocks{iB}).staAll.time, PPC.(Blocks{iB}).staAll.avg(:,:)');
%         ha = title(cfg_sta.spikechannel); set(ha,'Interpreter','none');
        set(gca,'FontSize',10,'XLim',cfg_sta.timwin,'XTick',cfg_sta.timwin(1):0.25:cfg_sta.timwin(2));
        xlabel('time (s)'); grid on;
    end
    %% ppc etc
    cfg_ppc                     = [];
    cfg_ppc.method              = 'mtmconvol';
    cfg_ppc.foi                 = 1:2:100;
    cfg_ppc.t_ftimwin           = 5./cfg_ppc.foi; % cycles per frequency
    cfg_ppc.taper               = 'hanning';
    cfg_ppc.spikechannel        = spk_chan;
    cfg_ppc.channel             = lfp_chan;
    cfg_ppc.rejectsaturation    = 'no'; % rejects saturation points which will cause PPC issues. 
    
    PPC.(Blocks{iB}).stsConvol  = ft_spiketriggeredspectrum(cfg_ppc, data_this_trl); % note, use raw or interpolated version
    
%     % plot
%     if isfield(cfg, 'plot')
%         figure(h)
%         subplot(2,2,2)
%         plot(PPC.(Blocks{iB}).stsConvol.freq,nanmean(sq(abs(PPC.(Blocks{iB}).stsConvol.fourierspctrm{1}))))
%     end
    %%
    cfg_ppc_stat                    = [];
    cfg_ppc_stat.method             = 'ppc0'; % compute the Pairwise Phase Consistency
    cfg_ppc_stat.spikechannel       = spk_chan;
    cfg_ppc_stat.channel            = lfp_chan;
    %cfg.dojack                = 1;
    cfg_ppc_stat.avgoverchan        = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
    cfg_ppc_stat.timwin             = 'all'; % compute over all available spikes in the window
    %cfg.latency               = [-2.5 0]; % sustained visual stimulation period
    PPC.(Blocks{iB}).statSts   = ft_spiketriggeredspectrum_stat(cfg_ppc_stat,PPC.(Blocks{iB}).stsConvol);
    
    % plot the results
    if isfield(cfg, 'plot')
        
        figure(iB)        
        subplot(222)
        plot(PPC.(Blocks{iB}).statSts.freq,PPC.(Blocks{iB}).statSts.ppc0')
        set(0,'DefaultTextInterpreter','none');
        set(gca,'FontSize',12);
        xlabel('frequency')
        ylabel('PPC')
        title(cfg_ppc.spikechannel);
    end
    obs_freq = PPC.(Blocks{iB}).statSts.freq;
    obs_ppc = PPC.(Blocks{iB}).statSts.ppc0';
    
    
    
    %%
    nShuf = cfg.shuffle;
    iChan = length(data_i.label);
%     data_i_shuff = data_this_trl; 
    
%     data_i.label{iChan} = 'temp_shuf';
    %%
    t_idx = strfind(data_this_trl.label, spk_chan);
    spk_idx = find(not(cellfun('isempty', t_idx)));
    shuf_ppc = zeros(nShuf,length(obs_freq));
    for iShuf = 1:nShuf
        fprintf('Shuffle %d...\n',iShuf);
        shuf_ppc(iShuf,:) =  Shuffle_PPC(cfg_ppc,cfg_ppc_stat, data_this_trl, spk_idx, lfp_chan, iChan);
    end
    
    id = strrep(spk_chan, '-', '_');
    
    z_ppc = (obs_ppc - nanmean(shuf_ppc,1)) / nanstd(shuf_ppc,1);
    
    %% plot
    if isfield(cfg, 'plot')
        figure(iB)        
        subplot(222)
        hold on;
        hx(1) = plot(obs_freq,obs_ppc,'k','LineWidth',2);
        plot(obs_freq,nanmean(shuf_ppc,1),'r');
        plot(obs_freq,nanmean(shuf_ppc,1)+nanstd(shuf_ppc,1),'r:');
        plot(obs_freq,nanmean(shuf_ppc,1)-nanstd(shuf_ppc,1),'r:');
        
        set(0,'DefaultTextInterpreter','none');
        legend(hx,{'observed','shuffled'},'Location','Northeast'); legend boxoff;
        set(gca,'FontSize',10);
        xlabel('frequency')
        ylabel('PPC')
        title([spk_chan '_' lfp_chan]);
        SetFigure([], gcf)
        
        %             mkdir(cfg.inter_dir, 'PPC')
            sess_id = strsplit(pwd, filesep);

        sess_id = strrep(sess_id{end}, '-', '_');

        saveas(gcf, [cfg.inter_dir filesep sess_id '_' id '_' lfp_chan '_' Blocks{iB}], 'fig');
        saveas(gcf, [cfg.inter_dir filesep  sess_id '_' id '_' lfp_chan '_' Blocks{iB}], 'png');

    end
    %% collect variables for export
    PPC.lfp_label = lfp_chan;
    PPC.spk_label = strrep(spk_chan, '-', '_');
    PPC.cfg_ppc = cfg_ppc;
    PPC.cfg_ppc_stat = cfg_ppc_stat; 
    PPC.cfg_sta = cfg_sta; 
    PPC.(Blocks{iB}).obs_freq = obs_freq;
    PPC.(Blocks{iB}).obs_ppc = obs_ppc;
    PPC.(Blocks{iB}).shuf_ppc = shuf_ppc;
    PPC.(Blocks{iB}).z_ppc = z_ppc;
    
    save([cfg.inter_dir  filesep 'PPC' sess_id '_' id '_' lfp_chan '.mat'], 'PPC', '-v7.3')
end


% if ~isempty(S_list)
%     [~,dir_id] = fileparts(pwd);
%     dir_id = strrep(dir_id, '-', '_');
%     dir_id = strrep(dir_id, ' ', '_');
%     save(['PPC_MS' dir_id], 'PPC',  '-v7.3');
% end

end
%%
function shuf_ppc = Shuffle_PPC(cfg_ppc,cfg_ppc_stat , data_i, spk_idx, lfp_chan, iChan)
% shuffle once
for iT = 1:length(data_i.trial) % shuffle each trial separately
    orig_data = data_i.trial{iT}(spk_idx,:);
    data_i.trial{iT}(iChan,:) = orig_data(randperm(length(orig_data)));
    
end

%% ppc etc
%     cfg_ppc                     = [];
%     cfg_ppc.method              = 'mtmconvol';
%     cfg_ppc.foi                 = 1:1:100;
%     cfg_ppc.t_ftimwin           = 5./cfg_ppc.foi; % cycles per frequency
%     cfg_ppc.taper               = 'hanning';
%     cfg_ppc.spikechannel        = spk_chan;
%     cfg_ppc.channel             = lfp_chan;
%     cfg_ppc.rejectsaturation    = 'no'; % rejects saturation points which will cause PPC issues. 
    
stsConvol     = ft_spiketriggeredspectrum(cfg_ppc , data_i); % note, use raw or interpolated version

% plot
%plot(stsConvol.freq,nanmean(sq(abs(stsConvol.fourierspctrm{1}))))

%%
% cfg_ppc_stat               = [];
% cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
% cfg_ppc.spikechannel = 'temp_shuf';
% cfg_ppc.channel      = lfp_chan;
% %cfg.dojack = 1;
% cfg_ppc.avgoverchan   = 'unweighted'; % weight spike-LFP phases irrespective of LFP power
% cfg_ppc.timwin        = 'all'; % compute over all available spikes in the window
% %cfg.latency       = [-2.5 0]; % sustained visual stimulation period
statSts           = ft_spiketriggeredspectrum_stat(cfg_ppc_stat ,stsConvol);

shuf_ppc = statSts.ppc0';


end % of shuffles
