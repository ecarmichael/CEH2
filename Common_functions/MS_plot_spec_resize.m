function [cut_vals, remove_flag] = MS_plot_spec_resize(cfg_in, data_in)
%% MS_plot_spec_resize: plots a spectrogram of the raw LFP and the traces for all LFP channels for
%   each segment in ms_seg.
%
%
%
%    Inputs:
%     - cfg_in: user configuration parameters
%
%     - data_in: [struct] contains all the LFP data.  [Alternative] use a
%     'tsd' structure like the output from MS_LoadCSC (WIP)
%
%    Outputs:
%     - cut_offs: [N x 2] array of user defined start and stop cut_offs
%     with the GUI.
%
%
%   ToDo: make this work for tsd data
%
% EC 2020-02-20   initial version
%   - ONly works for ms_seg at the moment.
%
%% inialize

cut_vals = NaN(2,length(data_in.time)); % fill in the cutoff values.
remove_flag = zeros(1,length(data_in.time));

cfg_def = [];
cfg_def.segments = 1:length(data_in.time); % which segments to process
cfg_def.emg_range = [-0.001 0.001]; % default, should be based on real data.
cfg_def.resize = 1; % use the gui to select the data for resizing. 
cfg_def.save_key = 'return'; % which key to use for saving the selected cutoffs
% cfg_def.redo_key = 'downarrow'; % which key for redoing the current cutoffs.
cfg_def.remove_key = 'backspace'; % which key to flag the current segment for removal.
cfg_def.spec.win_s = 2^10; % spectrogram window size.
cfg_def.spec.onverlap = cfg_def.spec.win_s / 2; % overlap
cfg_def.spec.freq = 10.5:0.1:80; % frequency range for spectrogram.
cfg_def.spec.lfp_chan = 2; % which channel to use for the spectrogram.  
cfg_def.saveas = [];

cfg = ProcessConfig(cfg_def, cfg_in);

if isfield('data_in', 'type') && strcmp(data_in.type, 'tsd')
    data_type = 'tsd';
else
    data_type = 'other';
end

%% loop through segments and plot.  If cfg.resize, then run the cutting GUI.

switch data_type
    
    case 'tsd'
        error('I haven''t built this...')
        
    case 'other'
        
        for iBlock  = cfg.segments % loop through segments.
            this_tvec = [];
            this_tvec = data_in.NLX_csc{iBlock}.tvec-data_in.NLX_csc{iBlock}.tvec(1);
            
            [~,F,T,P] = spectrogram(data_in.NLX_csc{iBlock}.data(cfg.spec.lfp_chan,:), rectwin(cfg.spec.win_s), cfg.spec.onverlap,cfg.spec.freq, data_in.NLX_csc{iBlock}.cfg.hdr{1}.SamplingFrequency);
            
            figure(iBlock+200)
            ax_spec(1) = subplot(7,1,1:3);
            ax1 = imagesc(T,F,10*log10(P));
            set(ax1, 'AlphaData', ~isinf(10*log10(P)))
            set(gca,'FontSize',10, 'xtick', []);
            axis xy; ylabel('Frequency (Hz)');
            ax = gca;
            % ax.YTickLabels = {0 [], [], 60, [], [], 80}; % set
            set(gca, 'tickdir', 'out');
            xlim([T(1) T(end)])

            title([data_in.hypno_label{iBlock} ' block id: ' num2str(iBlock)]);
            
            ax_spec(2) = subplot(7,1,4:6);
            hold on
            for iChan = length(data_in.NLX_csc{iBlock}.label):-1:1
                if strfind(data_in.NLX_csc{iBlock}.label{iChan}, '/')
                    plot(this_tvec, ((data_in.NLX_csc{iBlock}.data(iChan,:)/max(data_in.NLX_csc{iBlock}.data(iChan,:)))*10)+iChan*15);
                    text(this_tvec(1), iChan*15, '0')
                    yticks(iChan) = iChan*15;
                    yticks_label{iChan} = data_in.NLX_csc{iBlock}.label{iChan};
                    % add a baseline for reference.
                    l_h = hline(iChan*15, 'k');
                    l_h.Color = [0.8 0.8 0.8]; l_h.LineStyle = '--';
                elseif strcmp(data_in.NLX_csc{iBlock}.label{iChan}, 'EMG')
                    ax_spec(3) = subplot(7,1,7);
                    
                    plot(this_tvec, data_in.NLX_csc{iBlock}.data(iChan,:));
                    %             text(this_tvec(1)-2, mean(ms_seg_resize.NLX_csc{iBlock}.data(iChan,:)), ms_seg_resize.NLX_csc{iBlock}.label{iChan})
                    
                    ylim(cfg.emg_range)
                else
                    plot(this_tvec, (data_in.NLX_csc{iBlock}.data(iChan,:)*10000)+iChan*15);
                    %             text(this_tvec(1)-1, mean((ms_seg_resize.NLX_csc{iBlock}.data(iChan,:)*10000)+iChan*15), ms_seg_resize.NLX_csc{iBlock}.label{iChan})
                    
                    yticks(iChan) = iChan*15;
                    yticks_label{iChan} = data_in.NLX_csc{iBlock}.label{iChan};
                end
            end
            
            ax_spec(3) = subplot(7,1,7);
            y_label = get(gca, 'yticklabel');
            y_label{ceil(end/2), :} = 'EMG';
            y_label{1} = num2str(ax_spec(3).YAxis.TickValues(1));
            y_label{end} = num2str(ax_spec(3).YAxis.TickValues(end));
            set(gca, 'yticklabel', y_label)
            xlim([this_tvec(1) this_tvec(end)])
            xlabel('Time (s)')
            
            
            ax_spec(2) = subplot(7,1,4:6);
            yticks_label = yticks_label(~cellfun('isempty',yticks_label));
            set(gca, 'ytick', yticks,'yticklabel',yticks_label,  'xtick', [])
            xlim([this_tvec(1) this_tvec(end)])
            
            linkaxes(ax_spec, 'x');
            SetFigure([], gcf);
%             Resize_figure
            
            
            if cfg.resize ==1 % get data points for resizing.  First point is the start, second click is the end. Don't zoom in when making the cutoffs.
                input('Paused for inspection. Press any key to select cutoffs')
                
                fprintf('\n<strong>To save selection press: %s, to flag this seg for removal press: %s, to redo press: any key.</strong>\n', cfg.save_key, cfg.remove_key)
                break_out =0;
                while break_out ~= 1
                     [cut_x,~] = ginput(2);
                    hold on
                    v_ax(1) = vline(cut_x(1), 'r');
                    t_ax(1) = text(cut_x(1), F(1), 'Cutoff start', 'color', 'r','FontSize',14 );
                    v_ax(2) = vline(cut_x(2), 'r');
                    t_ax(2) = text(cut_x(2), F(1), 'Cutoff end', 'color', 'r','FontSize',14 );

                    was_a_key = waitforbuttonpress;
                    key_hit = get(gcf, 'CurrentKey');
                    if strcmp(key_hit, cfg.save_key)
                        remove_flag(iBlock) = 0;  
                        break_out = 1; 
                        if cut_x(1) <= this_tvec(1) && cut_x(2) >= this_tvec(end)
                            fprintf('\n<strong>TMS_plot_spec_resize:</strong> Segment #%d unchanged\n', iBlock)
                        else
                            fprintf('\n<strong>TMS_plot_spec_resize:</strong> Segment #%d flagged for resizing\n', iBlock)
                            
                        end
                        
                        
                    elseif strcmp(key_hit, cfg.remove_key)
                       remove_flag(iBlock) = 1;  
                       break_out = 1; 
                       fprintf('\n<strong>TMS_plot_spec_resize:</strong> Segment #%d flagged for removal\n', iBlock)
                       
                    else
                        fprintf('\n Redo...\n')

                    end
                    delete(v_ax);
                    delete(t_ax);

                end
                    
                if cut_x(1) <= this_tvec(1)
                    cut_x(1) = NaN;
                else
                    % get the time in nlx values.
                    [~,idx]=min(abs(this_tvec-cut_x(1)));
                    cut_x(1) =data_in.NLX_csc{iBlock}.tvec(idx);
                end
                
                if cut_x(2) >= this_tvec(end)
                    cut_x(2) = NaN;
                else
                    [~,idx]=min(abs(this_tvec-cut_x(2)));
                    cut_x(2) =data_in.NLX_csc{iBlock}.tvec(idx);
                end
                cut_vals(:,iBlock) = cut_x;
                
            else
                pause(2)
            end
            
            if ~isempty(cfg.saveas)
                if ~exist('Seg_plots', 'dir')
                mkdir('Seg_plots')
                end
                
                saveas(gcf, [pwd filesep 'Seg_plots' filesep 'Seg_' data_in.time_labels{iRec} '_' data_in.hypno_label{iRec}], cfg.saveas)
            end
                close
            
        end
        
end

%% clean up and output    
    