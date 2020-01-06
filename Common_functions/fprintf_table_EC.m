function fprintf_table_EC(fid,mat_in,rownames, colnames)
%% hack function to print a table as you would with fprintf
%
%
% inputs: 
%       fid: id from fopen
%       mat_in: your unput matrix
%       rownames: should be a cell array
%       col names: cell array
% example
%
%
%
%
%
%

if nargin <2
    error('need an FID and a matrix')
elseif nargin ==2
    warning('no row or column names.  filling in with blanks')
    colnames = [];
    rownames = [];
end
    

%% loop through rows and print

fprintf(fid, '\t');

for jj = 1:size(mat_in,2)
    fprintf(fid, ['\t' colnames{jj}]);
end

fprintf(fid, '\n\n');

% swap if the input matrix is a cell array
if ~iscell(mat_in)    
    format_str = repmat('\t%.3f',1,length(mat_in));

    for ii = 1:size(mat_in,1)
        fprintf(fid, [rownames{ii} '\t' format_str '\n'], mat_in(ii,:));      % Write Rows
    end
else
    format_str = repmat('\t%.3s',1,length(mat_in));

    for ii = 1:size(mat_in,1)
        fprintf(fid, [rownames{ii} '\t' format_str '\n'], mat_in{ii,:});      % Write Rows
    end
end
