function ms_out = MS_Remove_trace(cfg_in, ms_in)
%% MS_Remove_trace: remove specified traces from the all the data fields.
%
%
%
%    Inputs:
%     - cfg_in: [struct] user configurations
%           - cfg.remove_idx = [1xn] array
%           - cfg.data_type = 'RawTraces' could also be 'BinaryTraces' or 'Filtertraces'
%
%     - ms_data_in: [struct] the miniscope data structure following
%     preprocessing
%
%
%
%    Outputs:
%     - ms_data_out: [struct] with the specified traces
%
%
%
% EC 2020-01-28   initial version
%
%
%% initialize

cfg_def = [];
cfg_def.data_type = 'RawTraces';
cfg_def.user_fields = {}; % can add in other fields here.
cfg_def.user_fields_idx = []; % corresponding  index
cfg_def.user_fields_dim = []; % corresponding dimension

cfg_def.remove_idx = [];


cfg = ProcessConfig(cfg_def, cfg_in);
%%

ms_out = ms_in;

if isempty(cfg.user_fields) % if the user field is empty it will find all the fields with the number of cells and remove the cells from them.
    
    known_cell_num = size(ms_in.RawTraces,2);
    fields = fieldnames(ms_in);
    for iF = 1:length(fields)
        if isnumeric(ms_in.(fields{iF}))
            field_size = size(ms_in.(fields{iF}));
            
            cell_idx = find(field_size == known_cell_num,1);
            if ~isempty(cell_idx)
                fprintf('Removing traces in %s...\n', fields{iF})
                if cell_idx == 1 && length(field_size) == 1
                    ms_out.(fields{iF})(cfg.remove_idx) = [];
                    
                elseif cell_idx == 1 &&  length(field_size) == 2
                    ms_out.(fields{iF})(cfg.remove_idx,:) = [];
                    
                elseif cell_idx == 1 &&  length(field_size) == 3
                    ms_out.(fields{iF})(cfg.remove_idx,:,:) = [];
                    
                elseif cell_idx == 2 &&  length(field_size) == 2
                    ms_out.(fields{iF})(:,cfg.remove_idx) = [];
                    
                elseif cell_idx == 2 &&  length(field_size) == 3
                    ms_out.(fields{iF})(:,cfg.remove_idx,:) = [];
                    
                elseif cell_idx == 3
                    ms_out.(fields{iF})(:,:,cfg.remove_idx) = [];
                    
                else
                    error('Dealing with more dimensions than I had planned for')
                end
            end
        end
    end
    
else % only use the user fields
    
    if ~isequal(length(cfg.user_fields), length(cfg.user_fields_idx), length(cfg.user_fields_dim))
        error('User fields require a matched index and dimension to use for cell removal');
    else
        for iF = 1:length(cfg.user_fields)
            fprintf('Removing cells from %s\n', cfg.user_fields{iF});
            
            if cfg.user_fields_dim ==1 && cfg.user_fields_idx ==1
                ms_out.(cfg.user_fields{iF})(cfg.user_fields_idx) = [];
                fprintf('Size of %s is now %d\n', cfg.user_fields{iF}, size(ms_out.(cfg.user_fields{iF})))
                
            elseif cfg.user_fields_dim ==2 && cfg.user_fields_idx ==1
                ms_out.(cfg.user_fields{iF})(cfg.user_fields_idx,:) = [];
                fprintf('Size of %s is now %d x %d\n', cfg.user_fields{iF}, size(ms_out.(cfg.user_fields{iF}),1),size(ms_out.(cfg.user_fields{iF}),2))
                
            elseif cfg.user_fields_dim ==2 && cfg.user_fields_idx ==2
                ms_out.(cfg.user_fields{iF})(:,cfg.user_fields_idx) = [];
                fprintf('Size of %s is now %d x %d\n', cfg.user_fields{iF}, size(ms_out.(cfg.user_fields{iF}),1),size(ms_out.(cfg.user_fields{iF}),2))
                
            elseif cfg.user_fields_dim ==3 && cfg.user_fields_idx ==1
                ms_out.(cfg.user_fields{iF})(cfg.user_fields_idx,:,:) = [];
                fprintf('Size of %s is now %d x %d x %d\n', cfg.user_fields{iF}, size(ms_out.(cfg.user_fields{iF}),1),size(ms_out.(cfg.user_fields{iF}),2),size(ms_out.(cfg.user_fields{iF}),3))
                
            elseif cfg.user_fields_dim ==3 && cfg.user_fields_idx ==2
                ms_out.(cfg.user_fields{iF})(:,cfg.user_fields_idx,:) = [];
                fprintf('Size of %s is now %d x %d x %d\n', cfg.user_fields{iF}, size(ms_out.(cfg.user_fields{iF}),1),size(ms_out.(cfg.user_fields{iF}),2),size(ms_out.(cfg.user_fields{iF}),3))
                
            elseif cfg.user_fields_dim ==3 && cfg.user_fields_idx ==3
                ms_out.(cfg.user_fields{iF})(:,:,cfg.user_fields_idx) = [];
                fprintf('Size of %s is now %d x %d x %d\n', cfg.user_fields{iF}, size(ms_out.(cfg.user_fields{iF}),1),size(ms_out.(cfg.user_fields{iF}),2),size(ms_out.(cfg.user_fields{iF}),3))
            else
                error('Dealing with more dimensions than I had planned for...')
            end
        end
    end
    
end

ms_out.numNeurons = size(ms_out.RawTraces,2); 


