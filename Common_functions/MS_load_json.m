function j_out = MS_load_json(fname)

if nargin < 1
    j_files = dir('*.json');

    fid = fopen([fileparts(j_files.folder) filesep j_files.name]);
else
    fid = fopen(fname);
end

    
raw = fread(fid,inf);
str = char(raw');
j_out = jsondecode(str);