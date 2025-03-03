function ms = MS_Ca_good_cells(ms, keep_idx)
%% function MS_Ca_good_cells:  uses the 'keep_idx' from MS_currate_cells and removes them from each field. Needs the 'keep_idx' field as either part of the struct or as an input



%% initialize

if nargin < 2
    if isfield(ms, 'keep_idx')
        keep_idx = ms.keep_idx;
    else
        error('No keep_idx specific as an input or as a field of ms')
    end
end


%% loop over fields and remove the 'bad' cells


f_names = fieldnames(ms);

for iF = 1:length(f_names)
    
    if strcmp(f_names{iF}, 'keep_idx')
        continue
    else
        [s1, s2, s3] = size(ms.(f_names{iF}));
        
        
        if s1 == length(keep_idx)
            
            ms.(f_names{iF})(~keep_idx)  =[];
            
        elseif s2 == length(keep_idx)
            
            ms.(f_names{iF})(:,~keep_idx)  =[];
            
        elseif s3 == length(keep_idx)
            
            ms.(f_names{iF})(:,:,~keep_idx)  =[];
            
        end
    end
end







end

