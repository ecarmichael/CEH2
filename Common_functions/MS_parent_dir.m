function p_dir = MS_parent_dir(dir_in)
%%MS_parent_dir: outputs the parent directory for the current or a
%%specified directory. 


if nargin < 1
    dir_in = cd; 
end

parts = strsplit(dir_in, filesep); 

p_dir = strjoin(parts(1:end-1), filesep) ; 
