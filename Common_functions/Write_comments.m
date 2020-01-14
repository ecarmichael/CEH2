function Write_comments(fname, author)
%% Write_comments:prints a template for the comment section of a new function. If a name is 
%  specified as an input it will print based on that.
%
%
%
%    Inputs: 
%     - fname: [string] name of the function
%
%
%
%    Outputs: 
%     - none. 
%
%
%
%
% EC 2020-01-14   initial version 
%
%
%% print the comments

if nargin ==0
    fname = '';
    author = 'EC'
elseif nargin ==1
    author = 'EC';
end
fprintf('\n%%%% %s:\n%%\n%%\n%%', fname) % initial lines
fprintf('\n%%    Inputs: \n%%     -\n') % input lines
fprintf('%%\n%%\n%%') % input lines

fprintf('\n%%    Outputs: \n%%     -\n') % input lines
fprintf('%%\n%%\n') % input lines
fprintf('%%\n%%\n') % input lines
fprintf('%% %s %s   initial version \n%%\n%%\n%%%%\n',author, datestr(date, 'yyyy-mm-dd')) % input lines

end