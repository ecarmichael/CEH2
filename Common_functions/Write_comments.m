function Write_comments(fname, author, cfg)
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
    author = 'EC';
    cfg_def = [];
elseif nargin ==1
    author = 'EC';
        cfg_def = [];
elseif nargin == 2
        cfg_def = [];
elseif nargin == 3
    cfg_def = 1; 
end
if cfg_def == 1
fprintf('function %s(cfg_in)', fname) % initial lines
else
    fprintf('function %s', fname) % initial lines
end
fprintf('\n%%%% %s:\n%%\n%%\n%%', fname) % initial lines
if cfg_def == 1
    fprintf('\n%%    Inputs: \n%%    - cfg [struct]   configuration see the defaults below. \n%%\n%%    - \n') % input lines
else
    fprintf('\n%%    Inputs: \n%%    -\n') % input lines
end
fprintf('%%\n%%\n%%') % input lines

fprintf('\n%%    Outputs: \n%%    -\n') % input lines
fprintf('%%\n%%\n') % input lines
fprintf('%%\n%%\n') % input lines
fprintf('%% %s %s   initial version \n%%\n%%\n%%\n',author, datestr(date, 'yyyy-mm-dd')) % input lines
fprintf('%%%% initialize\n')
if cfg_def
fprintf('cfg_def = [];\n\n\n\n');
fprintf('cfg = ProcessConfig(cfg_def, cfg_in);\n%%%%\n');
end