function aligned_struct = MS_align_data_uni(struct2align, struct2match);
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

% get the time vector for the align2match
if isfield(struct2match, 'time')
    match_time_name = 'time'; 
elseif isfield(struct2match, 'tvec')
        match_time_name = 'tvec'; 
else
    error('needs a field in the struct2match that has time. Should be ''time'' or ''tvec''')
end


if isfield(struct2align, 'time')
    time_name = 'time'; 
elseif isfield(struct2align, 'tvec')
        time_name = 'tvec'; 
else
    error('needs a field in the struct2align that has time. Should be ''time'' or ''tvec''')
end



    [~, u_index] = unique(struct2align.(time_name)); % used to make sure there are no dubplicates which will break the interp1.


for iF =1:length(fieldnames) % loop through thr fields that are the same length as the time vector in struct2align

    if sum(ismember( length(struct2align.(time_name)), size(aligned_struct.(fieldnames{iF}))))
        if strcmp(fieldnames{iF}, 'vidNum') || strcmp(fieldnames{iF}, 'frameNum') 
            continue
        end
       
        fprintf('<strong>%s</strong>: interp on field: <strong>%s</strong>...\n', mfilename, fieldnames{iF})
        aligned_struct.(fieldnames{iF}) = [];
        
        if strcmp(fieldnames{iF}, time_name)
            aligned_struct.(fieldnames{iF}) = struct2match.(match_time_name);
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
                aligned_struct.(fieldnames{iF})(:,iD) = interp1(struct2align.(time_name)(u_index),struct2align.(fieldnames{iF})(u_index,iD),struct2match.(match_time_name));
            elseif dim == 1
                aligned_struct.(fieldnames{iF})(iD,:) = interp1(struct2align.(time_name)(u_index),struct2align.(fieldnames{iF})(iD,u_index),struct2match.(match_time_name));
            end
            
%             if isnan(aligned_struct.(fieldnames{iF})(end))
%                 aligned_struct.(fieldnames{iF})(end) = aligned_struct.(fieldnames{iF})(end-1)+mode(diff(    
%             end
        end
    end
    
end
