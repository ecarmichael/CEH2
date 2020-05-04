function h = MS_plot_ca_nlx(cfg_in, ms_data_in, csc)
%% MS_plot_ca_nlx: makes a simple plot of the NLX and the calicum traces. Can be either segmented (recording blocks) or one continous signal.
%
%
%
%    Inputs:
%     - cfg [struct] contains configuration paramters [see below]
%     - ms_data_in
%           if continuous, this is the ms.RawTraces or ms.FiltTraces [time x cell]
%           if segmented, this is a cell array with each containing
%           RawTraces or Filttraces [nSeg cells]
%     - csc :
%           if continuous, this is one csc signal [time x data]
%           if segmented, this is a cell array with each containing [time x
%           data]
%
%
%
%    Outputs:
%     - h : plot handle.
%
%       ToDO:
%           - switch based on ms_in being cell or not
%           - use csc channel labels to pick a channel when more than one
%           exist.
%
%
% EC 2020-01-20   initial version
%
%
%% initialize

% if the ms_data in is continous or segmented (in cells) set the type.
if ~iscell(ms_data_in.RawTraces)
    type = 'continuous';
    n_cells = size(ms_data_in.RawTraces,2);
    n_seg = 1;
    these_csc = csc.label; % get all the csc channel labels
else
    type = 'segmented';
    n_cells = size(ms_data_in.RawTraces{1},2);
    n_seg = length(ms_data_in.RawTraces); % how many recording segments.
    these_csc = csc{1}.label; 
end

if (iscell(ms_data_in.RawTraces) && ~iscell(csc)) || (~iscell(ms_data_in.RawTraces) && iscell(csc))
    error('Data are not consistent. Both ms_data_in and csc must both be continuous or segment, not a mix')
end

% get number of Ca cells.



%% set up defaults.
cfg_def =[];
cfg_def.Ca_type = 'RawTraces'; % can be either 'RawTraces' or FiltTraces' (maybe others?)
cfg_def.Ca_chan = 1:floor(n_cells/10):n_cells; % get a subset of cells.
cfg_def.chan_to_plot = these_csc;  % how many channels are in the csc.  Determined above. 
cfg_def.plot_type = '2d'; % '2d' or '3d'
cfg_def.x_zoom = []; % where to zoom in on eht x_axis for each plot.
cfg_def.emg_range = [];
cfg_def.label = 'file_names'; 
cfg_def.saveas = []; % if empty don't save the images.  Can be '.fig' (matlab fig) or other known saveas format.  I like '.png'. 

cfg = ProcessConfig(cfg_def, cfg_in);


%% make the plots

c_ord = linspecer(length(cfg.Ca_chan)); % nice colours.
c_ord_lfp = linspecer(length(cfg.chan_to_plot));
switch type
    
    case 'segmented'
        
        for iRec = n_seg:-1:1
            h(iRec) = figure(iRec);
            timein = (csc{iRec}.tvec - csc{iRec}.tvec(1)); % just to fix the timing offset between them back to ebing relative to this segment.
            for iChan = 1:length(cfg.chan_to_plot)
                ax(iChan) =subplot(length(cfg.chan_to_plot)+2,1,iChan);
                
                plot(timein,csc{iRec}.data(iChan,:), 'color', c_ord_lfp(iChan,:));
                if strcmp(csc{iRec}.label{iChan}, 'EMG') && ~isempty(cfg.emg_range)
                    ylim(cfg.emg_range)
                end
                title(csc{iRec}.label{iChan})
                xlim([timein(1), timein(end)])
            end
            
            ax(length(cfg.chan_to_plot)+1) =subplot(length(cfg.chan_to_plot)+2,1,[length(cfg.chan_to_plot)+1 length(cfg.chan_to_plot)+2]);
            hold on
            time_in2 = ms_data_in.time{iRec} - ms_data_in.time{iRec}(1);
            switch cfg.plot_type
                % 2d
                case '2d'
                    for iC = 1:length(cfg.Ca_chan)
%                         if strcmp(cfg.Ca_type, 'BinaryTraces')
                            plot(time_in2*0.001, ms_data_in.(cfg.Ca_type){iRec}(:,iC)+iC, 'color', c_ord(iC,:))
%                             area(time_in2*0.001, max(ms_data_in.(cfg.Ca_type){iRec}(:,iC)+iC, iC), iC, 'edgecolor', 'none', 'facecolor', c_ord(iC,:))
%                         else
%                             plot(time_in2*0.001, ms_data_in.(cfg.Ca_type){iRec}(:,iC), 'color', c_ord(iC,:))
%                         end
                    end
                    % 3d
                case '3d'
                    for iC = 1:length(cfg.Ca_chan)
                        plot3(time_in2*0.001, repmat(iC,size(ms_data_in.(cfg.Ca_type){iRec},1),1), ms_data_in.(cfg.Ca_type){iRec}(:,iC)*3, 'color', c_ord(iC,:))
                    end
                    view([0 75])
                    set(gca, 'color', [.5 .5 .5])
            end
            if strcmp(cfg.Ca_type, 'BinaryTraces')
                set(gca, 'yticklabel',cfg.Ca_chan);
                ylabel('Cell ID');
            end
            xlim([time_in2(1)*0.001 time_in2(end)*0.001])
            xlabel('time (s)')
%             linkaxes(ax, 'x') % only meant for linking 2d axes. 
            link = linkprop(ax,{'XLim'}); % links x axes across 2d and 3d plots. 
            setappdata(gcf, 'StoreTheLink', link)
%             ax = [];
            if ~isempty(cfg.x_zoom)
            xlim(cfg.x_zoom)
            end
            title(ms_data_in.(cfg.label){iRec})
            if ~isempty(cfg.saveas)
                mkdir('Seg_plots')
                
                saveas(gcf, [pwd filesep 'Seg_plots' filesep 'Seg_' ms_data_in.time_labels{iRec} '_' ms_data_in.hypno_label{iRec} cfg.saveas])
            end
        end
        
    case 'continuous'
        warning('Plotting continuous Ca data does not make sense ATM since the CA is in blocks while the CSC is not')
        h = figure(1);
        ax(1) =subplot(2,1,1);
        timein = (csc.tvec - csc.tvec(1)); % just to fix the timing offset between them back to ebing relative to this segment.
        
        plot(timein,csc.data, '-b' );
        xlim([timein(1), timein(end)])
        
        ax(2) =subplot(2,1,2);
        hold on
        time_in2 = ms_data_in.time - ms_data_in.time(1);
        switch cfg.plot_type
            % 2d
            case '2d'
                for iC = 1:length(cfg.Ca_chan)
                    plot(time_in2*0.001, ms_data_in.(cfg.Ca_type)(:,iC), 'color', c_ord(iC,:))
                end
                % 3d
            case '3d'
                for iC = 1:length(cfg.Ca_chan)
                    plot3(time_in2*0.001, repmat(iC,size(ms_data_in.(cfg.Ca_type),1),1), ms_data_in.(cfg.Ca_type)(:,iC), 'color', c_ord(iC,:))
                end
                view([0 45])
        end
        
        xlim([time_in2(1)*0.001 time_in2(end)*0.001]);
        linkaxes(ax, 'x');
%         ax = [];
        xlim(cfg.x_zoom);
        
end % end swtich
end % end function