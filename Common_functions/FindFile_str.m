function f_id = FindFile_str(dir_name, str)
%% FindFile_str: 
%   finds all files in a directory with a specified string "str" in the filename
%
%
% input:
%       dir: [string] with the path to the directory you wish to search
%       str: [string] string to find
%
% outputs: 
%       f_id: [cell] contains the file names that contained the string.  
%
% found on matlab centeral: https://www.mathworks.com/matlabcentral/answers/101485-how-do-i-find-a-file-containing-a-particular-string-in-a-given-directory-in-matlab-7-13-r2011b

% modified into a function by EC 2018

filesAndFolders = dir(dir_name);     % Returns all the files and folders in the directory
filesInDir = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory                    
stringToBeFound = str;
numOfFiles = length(filesInDir);
ii=1;
f_id = {[]}; 
while(ii<=numOfFiles)
      filename = filesInDir(ii).name;                              % Store the name of the file
%       while(~feof(fid))                                           % Execute till EOF has been reached
found = regexp( filename ,regexptranslate('wildcard',stringToBeFound)); % better way to deal with wildcards in the middle of string
% found = strfind(filename,stringToBeFound);         % Search for the stringToBeFound in contentOfFile
          if ~isempty(found)
              foundString = strcat('Found in file------', filename);
              disp(foundString);
              f_id{end+1} = filename;
          end   
            
      ii = ii+1;
end
f_id = f_id(~cellfun('isempty',f_id));
end