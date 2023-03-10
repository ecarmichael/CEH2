function block_idx = MS_get_waze_blocks(evt)

%% MS_get_waze_blocks:
%
%
%
%    Inputs: 
%    -
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2022-09-20   initial version 
%
%
%
%% initialize
start_idx = find(contains(evt.label, 'Starting Recording'));
end_idx = find(contains(evt.label, 'Stopping Recording'));




if length(evt.t{start_idx}) ~= 4
    rec_len = evt.t{end_idx} - evt.t{start_idx};
    [~, idx] = sort(rec_len, 'descend'); 
    keep_idx = sort(idx(1:4)); 
    
%     error('too many start and stop times. should be 4.  Figure this out')
else
    keep_idx = 1:4; 
end

block_idx = [];
if length(evt.t{start_idx}) > 1
    disp('splitting into blocks')
    Blocks = {'Pre_sleep', 'W_maze', 'OF', 'Post_sleep'};
    for iB = 1:length(Blocks)
        if strcmp(Blocks{iB}, 'Pre_sleep')  || strcmp(Blocks{iB}, 'Post_sleep') % trim sleep blocks. 
            block_idx.(Blocks{iB}) = [evt.t{start_idx}(keep_idx(iB)) evt.t{start_idx}(keep_idx(iB))+120*60];
        else
        block_idx.(Blocks{iB}) = [evt.t{start_idx}(keep_idx(iB)) evt.t{end_idx}(keep_idx(iB))];
        end
        fprintf('Block: %s duration = %.0fmins\n', Blocks{iB}, (block_idx.(Blocks{iB})(2) -block_idx.(Blocks{iB})(1))/60); 
    end
end
