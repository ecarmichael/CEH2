function ms_seg_out = MS_resize_segments(cfg_in, ms_seg_in)
%% MS_resize_Segments: resizes the data segments in the ms_seg structure given cfg.cutoff values.
%       If cfg.cutoffs are in the Neuralynx time format (cfg.tvec_to_use =
%       'NLX_csc') then the NLX_evt (events) field will be used to find the
%       nearest index in the miniscope data and the NLX_csc structure will
%       be restricted using 'restrict.m'.  If the cutoffs are in 'miniscope
%       time' (cfg.tvec_to_use = 'time') then the exact indices will be
%       used for restricting miniscope data fields and the NLX_csc fields
%       will be determined using the nearest values. 
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

switch cfg.tvec_to_use
    
    case 'NLX_csc'
        % find the index in the NLX_evt.t that is the closest to the
        % nlx_csc cutoff times.  Use the actual time for nlx_csc which uses
        % restrict to save the time of cycling through fields. 
        for iC = length(cfg.cutoffs):-1:1
            if isnan(cfg.cutoffs(1,iC))
                tstart_idx.time(iC) = 1;
                tstart_idx.nlx(iC) = ms_seg_in.NLX_csc{iC}.tvec(1);
            else
                tstart_idx.time(iC) = nearest_idx3(cfg.cutoffs(1,iC), ms_seg_in.NLX_evt{iC}.t{end});
                tstart_idx.nlx(iC) = ms_seg_in.NLX_csc{iC}.tvec(nearest_idx3(cfg.cutoffs(1,iC), ms_seg_in.NLX_csc{iC}.tvec));
            end
            
            if isnan(cfg.cutoffs(2,iC))
                tend_idx.time(iC) = length(ms_seg_in.NLX_evt{iC}.t{end});
                tend_idx.nlx(iC) =ms_seg_in.NLX_csc{iC}.tvec(end);

            else
                tend_idx.time(iC) = nearest_idx3(cfg.cutoffs(2,iC), ms_seg_in.NLX_evt{iC}.t{end});
                tend_idx.nlx(iC) = ms_seg_in.NLX_csc{iC}.tvec(nearest_idx3(cfg.cutoffs(2,iC), ms_seg_in.NLX_csc{iC}.tvec));
            end
            
        end        
        
    case 'time'
         for iC = length(cfg.cutoffs):-1:1
            if isnan(cfg.cutoffs(1,iC))
                tstart_idx.time(iC) = 1;
            else
                tstart_idx.time(iC) = cfg.cutoffs(1,iC);
            end
            
            if isnan(cfg.cutoffs(2,iC))
                tend_idx.time(iC) = length(ms_seg_in.time{iC});

            else
                tend_idx.time(iC) = cfg.cutoffs(2,iC);
            end
            
        end      
        
        
%         error('Haven''t built this yet. Probably just need to copy above but look up nlx_evt. EC 2020-02-18')

end


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
            
            for iSeg = length(ms_seg_in.(cfg.tvec_to_use)):-1:1
                if sum(isnan(cfg.cutoffs(:,iSeg))) ==2
                    ms_seg_out.(fields{iF}){iSeg} = ms_seg_in.(fields{iF}){iSeg};
                else
                    if strcmp(fields{iF}, 'NLX_csc') || strcmp(fields{iF}, 'NLX_evt')
                       ms_seg_out.(fields{iF}){iSeg} = restrict(ms_seg_out.(fields{iF}){iSeg},tstart_idx.nlx(iSeg), tend_idx.nlx(iSeg));

                    else
                        ms_seg_out.(fields{iF}){iSeg} = [];
                        ms_seg_out.(fields{iF}){iSeg} = ms_seg_in.(fields{iF}){iSeg}(tstart_idx.time(iSeg): tend_idx.time(iSeg));

                    end
                end
            end
        end
    end
end
                

%% clean up
ms_seg_out.resize.cfg = cfg;
ms_seg_out.resize.cfg.date = date; 


end


