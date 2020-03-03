function ExpKeys = MS_Load_Keys(file_name)
%% MS_Load_keys: loading function for the 'Keys.m' files.  Keys files are 
% used to track all experimental variables like subject, tragets, behaviour...
% 
%
%
%   Inputs:
%       - none: will look for anything ending in '*Keys.m'
%       - file_name: [string] file to run as a string (eg:
%       'M537day0base1_Keys.m', or
%       '/Users/jericcarmichael/Documents/Williams_Lab/M537day0base1_Keys.m').
%        
%   Outputs:
%       - none.  Instead it will run the *Keys.m file resulting in the Keys
%       variable. 
%
%
%   EC 2019-12-30: initial version based on ExpKeys/LoadExpKeys/FindFile by
%   MvdM and AD Redish
%
%% extract inputs and run or 

if nargin ==1 % if specified, just run it. 
    run(file_name); 
    fprintf(['\nKeys file: ' file_name ' loaded.\n']); 
    
else
    file_name = dir('*Keys.m');
    if isempty(file_name)
        error('No *Keys.m files found here')
    end
        run(file_name.name); 
    fprintf(['\nKeys file: ' file_name.name ' loaded.\n']); 
    

end
ExpKeys = Keys; 
clear file_name
end