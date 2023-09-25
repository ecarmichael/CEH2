function aligned_struct = MS_align_data(struct2align, struct2match);
%% MS_align_data: aligns one miniscope struct (struct2align) and a template (struct2match).
%       MS_align_data will cycle through all the continuous sample fields
%       in struct2align and interpolate them to match the length of the
%       time vector in struct2match.  Key fields may be behav.position, behav.speed,...
%
%
%
%    Inputs:
%    - struct2align: [struct] contains either miniscope ('ms') or
%    behavioural ('behav') structures from msNormCorre or msExtractBehavoir
%
%
%
%    Outputs:
%    - algined_struct: [struct] interpolated version of struct2align
%
%
%
%
% EC 2020-07-14   initial version
%
%% cycle through fields in Struct2algin and interp to match struct2match.time
aligned_struct = struct2align; % clone the struct

fieldnames = fields(struct2align);

[~, u_index] = unique(struct2align.time); % used to make sure there are no dubplicates which will break the interp1.


for iF =1:length(fieldnames) % loop through thr fields that are the same length as the time vector in struct2align

    if length(aligned_struct.(fieldnames{iF})) == length(struct2align.time)
        if strcmp(fieldnames{iF}, 'vidNum') || strcmp(fieldnames{iF}, 'frameNum') 
            continue
        end
       
        fprintf('<strong>%s</strong>: interp on field: <strong>%s</strong>...\n', mfilename, fieldnames{iF})
        aligned_struct.(fieldnames{iF}) = [];
        
        if strcmp(fieldnames{iF}, 'time') 
            aligned_struct.(fieldnames{iF}) = struct2match.time;
            continue
        end
        %         if size(aligned_struct.(fieldnames{iF}),2) == 1
        [m, n] = size(struct2align.(fieldnames{iF}));
        if m > n
            dim = 2;
        else
            dim = 1;
        end
        for iD = 1:size(struct2align.(fieldnames{iF}),dim)
            if dim == 2
                aligned_struct.(fieldnames{iF})(:,iD) = interp1(struct2align.time(u_index),struct2align.(fieldnames{iF})(u_index,iD),struct2match.time);
            elseif dim == 1
                aligned_struct.(fieldnames{iF})(iD,:) = interp1(struct2align.time(u_index),struct2align.(fieldnames{iF})(iD,u_index),struct2match.time);
            end
            
%             if isnan(aligned_struct.(fieldnames{iF})(end))
%                 aligned_struct.(fieldnames{iF})(end) = aligned_struct.(fieldnames{iF})(end-1)+mode(diff(    
%             end
        end
    end
    
end
