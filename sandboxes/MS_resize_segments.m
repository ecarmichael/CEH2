function ms_seg_out = MS_resize_Segments(cfg_in, ms_seg_in)
%% MS_resize_Segments: resizes the data segments in the ms_seg structure given cfg.cutoff values.
%
%
%
%    Inputs:
%     - cfg_in: [struct] contains user configations.
%           - cfg.cutoffs = [2 x nSegments] array.  Contains the start and
%               stop cut offs (default: based on the NLX_csc.tvec time values). If NaN vals, these will output the full
%               segment.
%           - cfg.tvec_to_use: [string] what channel do you want to base
%           the cutoffs on? typically ms_seg.NLX_csc.tvec, but could also be
%           ms_seg.time
%
%     - ms_seg_in: [struct] contains the miniscope data which has been
%           segmented (MS_segment_ms_sandbox).
%
%
%    Outputs:
%     - ms_seg_out: [struct] same format as the ms_seg input but with the new
%           restricted segments.
%
%
%   ToDo: make this work with just the ms time values so that it can be
%   indpendent of the NLX_csc data.
%
% EC 2020-02-14   initial version
%
%
%% initialize
cfg_def = [];
cfg_def.tvec_to_use = 'time'; % what subfield to use as the tvec/time scale (default is ms_seg.time using miniscope timestamps).
cfg_def.cutoffs = NaN(length(ms_seg_in.(cfg_def.tvec_to_use)), 2); % this would result in nothing being restricted.

cfg = ProcessConfig(cfg_def, cfg_in);

% if strcmp(cfg.tvec_to_use, 'NLX

fprintf('\n<strong>MS_resize_Segments</strong>: resizing segments using the %s time scale\n', cfg.tvec_to_use);

%% loop through segments and resize all the data in the ms_seg_in
ms_seg_out = ms_seg_in;




known_cell_num = size(ms_seg_in.RawTraces,1); % should always be the correct number of cells for the number of segments.
fields = fieldnames(ms_seg_in);
for iF = 1:length(fields)
    field_size = size(ms_seg_in.(fields{iF}));
    
    cell_idx = find(field_size == known_cell_num,1);
    
    
    if ~isempty(cell_idx)
        if ischar(ms_seg_in.(fields{iF}){1})
            continue
        else
            fprintf('Resizing traces in <strong>''%s''</strong>...\n', fields{iF})
            %             if cell_idx == 1 && length(field_size) == 1
            %                 ms_out.(fields{iF})(cfg.remove_idx) = [];
            %
            %             elseif cell_idx == 1 &&  length(field_size) == 2
            %                 ms_out.(fields{iF})(cfg.remove_idx,:) = [];
            %
            %             elseif cell_idx == 1 &&  length(field_size) == 3
            %                 ms_out.(fields{iF})(cfg.remove_idx,:,:) = [];
            %
            %             elseif cell_idx == 2 &&  length(field_size) == 2
            
            for iSeg = length(ms_seg_in.(cfg.tvec_to_use)):-1:1
                if sum(isnan(cfg.cutoffs(:,iSeg))) ==2
                    ms_seg_out.(fields{iF}){iSeg} = ms_seg_in.(fields{iF}){iSeg};
                else
                    if strcmp(cfg.tvec_to_use, 'NLX_csc') || strcmp(cfg.tvec_to_use, 'NLX_evt')
                        
                        if strcmp(fields{iF}, 'NLX_csc') || strcmp(fields{iF}, 'NLX_evt')
                        ms_seg_out.(fields{iF}){iSeg} = restrict(ms_seg_in.(fields{iF}){iSeg}, cfg.cutoffs(1,iSeg), cfg.cutoffs(2,iSeg));
                        
                        
                        fprintf('    Resizing segment # %d between %0.1fs - %0.1fs\n', iSeg, ms_seg_in.(fields{iF}){iSeg}.tvec(1) - ms_seg_in.(fields{iF}){iSeg}.tvec(idx(1)), ms_seg_in.(fields{iF}){iSeg}.tvec(end) - ms_seg_in.(fields{iF}){iSeg}.tvec(idx(2)))
                        else
                            [~, idx] =intersect(ms_seg_in.(cfg.tvec_to_use){iSeg}.tvec,cfg.cutoffs(:,iSeg));
                            this_tvec = ms_seg_in.(cfg.tvec_to_use){iSeg}.tvec - ms_seg_in.(cfg.tvec_to_use){iSeg}.tvec(1);
                            
                            
                            
                            ms_seg_out.(fields{iF}){iSeg} = restrict(ms_seg_in.(fields{iF}){iSeg}, cfg.cutoffs(1,iSeg), cfg.cutoffs(2,iSeg));

                            
                            
                        end
                        else
                        [~, idx] =intersect(ms_seg_in.(cfg.tvec_to_use){iSeg},cfg.cutoffs(:,iSeg));
                        ms_seg_out.(fields{iF}){iSeg} = ms_seg_in.(fields{iF}){iSeg}(idx(1):idx(2));
                        fprintf('    Resizing segment # %d between %0.1fs - %0.1fs\n', iSeg, ms_seg_in.(fields{iF}){iSeg}(1) - ms_seg_in.(fields{iF}){iSeg}(idx(1)), ms_seg_in.(fields{iF}){iSeg}(end) - ms_seg_in.(fields{iF}){iSeg}(idx(2)))
                        
                    end
                end
            end
            
            
            %             elseif cell_idx == 2 &&  length(field_size) == 3
            %                 ms_out.(fields{iF})(:,cfg.remove_idx,:) = [];
            %
            %             elseif cell_idx == 3
            %                 ms_out.(fields{iF})(:,:,cfg.remove_idx) = [];
            %
            %             else
            %                 error('Dealing with more dimensions than I had planned for')
            %             end
        end
    end
    
end

%% clean up



end


