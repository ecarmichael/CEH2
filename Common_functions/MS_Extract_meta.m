function META = MS_Extract_meta(fname, save_dir)
%% MS_Extract_meta: reads .json file from the 
%
%
%
%    Inputs: 
%    - fname [string or path] OPTIONAL the name of the file to load. If
%    empty will look for a 'metaData.json' file in cd. 
%
%    save_dir:  [string or path]  OPTIONAL location to save the META.m
%    file. If empty it will save to the cd. 
%
%    Outputs: 
%    - META  [struct]  contains all the metadata from the .json file and
%    any other defaults infered from the .json. 
%
%
%
%
% EC 2021-05-28   initial version 
%
%
%
%% initialize

if nargin <1
    fname = 'metaData.json'; 
end

if nargin < 2
    save_dir = cd; 
end

% check for the file
if ~exist(fname, 'file')
    error('''%s'' Not found', fname)
end

%% open and collect fields from .json file. 
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
META = jsondecode(str);

META.date = datestr([num2str(META.year) '-' datestr(META.month,'mm') '-' datestr(META.day, 'dd')], 'yyyy_mm_dd');
META.time = [num2str(META.hour,'%02.f') ':' num2str(META.minute,'%02.f') ':' num2str(META.second,'%02.f') '.' num2str(META.msec,'%03.f')];
META = rmfield(META, {'hour', 'minute', 'second', 'year', 'month', 'day', 'msec'});

%% add some other fields 

if contains(META.animalName, 'ck2')
    META.geno = 'ck2';
elseif contains(META.animalName, 'nt')
    META.geno = 'nt';
end

%% convert META in an exicutible script

% open the file to write
fid = fopen([save_dir filesep 'Meta.m'], 'w');

% fill in the consistent information
fprintf(fid, ['%% This META.m was generated using MS_Extract_meta.m on ' num2str(date) ';\n']);

fprintf(fid, ['META.version = ' num2str(1) ';\n']);


these_fields = fieldnames(META);
for iF = 1:length(these_fields)
    
    
    
if ~iscell(META.(these_fields{iF}))
    if ischar(META.(these_fields{iF}))
        fprintf(fid, ['META.' these_fields{iF} ' = ''' META.(these_fields{iF}) ''' ;\n']);
    else
        fprintf(fid, ['META.' these_fields{iF} ' = ' num2str(META.(these_fields{iF})) ';\n']);
    end
end
        
    
    
    
end

fprintf(fid, '\n%%Notes\n');
fprintf(fid, 'META.notes = '''';\n');

fclose(fid);
disp('META written')
%%