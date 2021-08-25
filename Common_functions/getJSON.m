function out = getJSON(fname)
%% Load and parse a .json file.  from MATLAB answers: https://www.mathworks.com/matlabcentral/answers/474980-extract-info-from-json-file-by-matlab


str = fileread(fname); % dedicated for reading files as text

out = jsondecode(str);
end