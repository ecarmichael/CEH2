function params = OE_load_params(params_dir)

if nargin < 1
    params_dir = cd; 
end

% find the params file
fname = dir([params_dir filesep 'params.py']); 

% read the file
fid = fopen([fname.folder filesep fname.name]); 
Filedata = textscan(fid, '%s%s', 'Delimiter', '=',  'ReturnOnError', false);
fclose(fid);

names = Filedata{1};  val = Filedata{2}; 

params = []; 
for ii = 1:length(names)
    if ~isnan(str2double(val{ii}))
        params.(names{ii}(1:end-1)) = str2double(val{ii});
    else
        params.(names{ii}(1:end-1)) = val{ii};
    end
end

params.dat_path = params.dat_path(3:end-2); 