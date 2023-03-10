function Minian_Merge(dir_list, merge_dir, scope_name)
%% Minian_Merge: collects all avi files in dir1_list and puts them in merge_dir along with the timestamps. Renames dir2 files to follow dir1 numerically. 
%
%
%
%    Inputs: 
%    - dir_list: [path] list of dirs to copy files from.  If this is 
%
%    - scope_name: [string]  name of the scope folder. 
%
%    -merge_dir: [path]    where to store the new files
%
%    Outputs: 
%    none
%
%
%
%
% EC 2023-02-05   initial version 
%
%  TO DO: make it check if dir1 is a list, if then 
%
%% initialize

if nargin <3
    scope_name = 'My_V4_Miniscope';
end

%% list all the dirs in the list and find the avi files
dir2merge = []; avi_list = []; ts_list = [];
for iD = 1:length(dir_list)
    dir2merge{iD} = [dir_list{iD} filesep scope_name];
%     cd([dir_list{iD} filesep scope_name])
    
    avi_list{iD} = dir([dir_list{iD} filesep scope_name filesep '*.avi']);
    fprintf('<strong>%.0f</strong> .avi files found in <strong> %s</strong>\n', length(avi_list{iD}), dir_list{iD}); 
    ts_list{iD} = dir([dir_list{iD} filesep scope_name filesep 'timeStamps.csv']);

end

%% move the files to ther merge dir and append 100s to maintain order
parts = strsplit(dir_list{1}, filesep);
sub = parts{contains(parts, 'dSub')}; 
sess =  parts{contains(parts, '2023')}; 

merge_dir_name = [merge_dir filesep sub '_' sess '_' 'merged'];
mkdir(merge_dir_name);

for iD = 1:length(avi_list)
   for ii = 1:length(avi_list{iD})
       
       vid_num = strsplit(avi_list{iD}(ii).name, '.');
       if length(vid_num{1}) > 1
           copyfile([avi_list{iD}(ii).folder filesep avi_list{iD}(ii).name],[merge_dir_name filesep num2str(iD*100) avi_list{iD}(ii).name]);
       else
           copyfile([avi_list{iD}(ii).folder filesep avi_list{iD}(ii).name],[merge_dir_name filesep num2str(iD*100) '0' avi_list{iD}(ii).name])
       end
   end
end


%% combine the TS files
TS = []; 
all_TS = [];
for iT = 1:length(ts_list)
    TS{iT} = readtable([ts_list{iT}.folder filesep ts_list{iT}.name]);
    all_TS = [all_TS;TS{iT}];
end

writetable(all_TS, [merge_dir_name filesep 'timeStamps.csv'])
