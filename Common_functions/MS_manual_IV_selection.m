function [keep_idx] = MS_manual_IV_selection(csc, iv_in, win)

% init
if nargin < 3
    win = .5; % .5 second window on either end. 
end

keep_idx = zeros(length(iv_in.tstart), 1); 

fprintf('Press <strong>[space bar]</strong> to keep; press any other to skip\n')

fig = figure(19191919); 

clf
cfg.target = csc.label{1}; 
PlotTSDfromIV(cfg, iv_in, csc)
% plot(csc.tvec, csc.data(1,:), 'k')
% hold on
% plot(csc.tvec, csc.data(1,:), 'k')
fprintf('\n\n')
last_s = []; 

for ii = 1:length(iv_in.tstart)

    s = sprintf('%04d / %d', ii, length(iv_in.tstart)); 

fprintf([repmat('\b', 1, length(last_s)) s])
last_s = s; 

xlim([iv_in.tstart(ii) - win  iv_in.tend(ii)+win])

    w = waitforbuttonpress; 
    
    % If a key was pressed, check if it was the space bar
    if w == 1 && strcmp(get(fig, 'CurrentCharacter'), ' ')
        keep_idx(ii) = true;         
    elseif  w == 1 && strcmp(get(fig, 'CurrentCharacter'), 'escape')
        return
    end

end