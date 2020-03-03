function ms_out = MS_de_cell(ms_in);
%% MS_de_cell:
%
%
%
%    Inputs: 
%     -
%
%
%
%    Outputs: 
%     -
%
%
%
%
% EC 2020-02-28   initial version 
%
%
%%  Loop through fields and if they don't need to be a cell convert them back to an array. 
ms_out = ms_in; 

known_cell_num = size(ms_in.RawTraces,1); % should always be the correct number of cells for the number of segments.
fields = fieldnames(ms_in);
for iF = 1:length(fields)
    field_size = size(ms_in.(fields{iF}));
    
    cell_idx = find(field_size == known_cell_num,1);
    
    
    if ~isempty(cell_idx) && iscell(ms_in.(fields{iF}))
        
            ms_out.(fields{iF}) = ms_in.(fields{iF}){1};
            
            
            fprintf('Converting Cells to array in %s\n', fields{iF});
    end
end

fprintf('<strong>%s</strong>: Converted Cells in data fields to array.\n', mfilename);

%% 