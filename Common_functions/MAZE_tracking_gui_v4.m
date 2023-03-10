function pos = MAZE_tracking_gui_v4(pos)
%% MAZE_tracking_gui: loads the position and event data (NLX format)
%
%
%
%    Inputs:
%    - pos [struct] position data in the TSD format
%
%    Outputs:
%    - maze_evts: [struct] contains outputs from the MAZE_gui
%
%
%
%
% EC 2022-09-15   initial version
%
%
%

%% run the MAZE GUI to get User defined experimental time points.

Maze_GUI_NLX_sandbox

% waitfor(0, 'Userdata', 'GUI_states');
%% convert to trials

% box_in = GUI_state_times(contains(GUI_states,  'Box_in'));
% states = unique(GUI_states);
% 
% maze = [];
% maze.trials.times = []; maze.trials.types = [];
% for ii = length(states):-1:1
%     
%     if contains(states{ii}, {'Box_in', 'F', 'C'})
%         this_evt = find(contains(GUI_states, states{ii}));
%         if this_evt(end) == length(GUI_state_times)
%             this_evt = this_evt(1:end-1);
%         end
%         
%         maze.events.(states{ii}) = [GUI_state_times(this_evt);  GUI_state_times(this_evt+1)]';
%         
%         if strcmp(states{ii}, 'C_R') || strcmp(states{ii}, 'C_L')
%             maze.trials.times = [maze.trials.times; maze.events.(states{ii})];
%             for iT = 1:length(this_evt)
%                 if  strcmp(GUI_states{this_evt(iT)}, 'C_R') && strcmp(GUI_states{this_evt(iT)+1}, 'R_Rew')
%                     maze.trials.types{end+1} = 'CR correct';
%                 elseif strcmp(GUI_states{this_evt(iT)}, 'C_R') && strcmp(GUI_states{this_evt(iT)+1}, 'L_Rew')
%                     maze.trials.types{end+1} = 'CR error';
%                 elseif  strcmp(GUI_states{this_evt(iT)}, 'C_L') && strcmp(GUI_states{this_evt(iT)+1}, 'L_Rew')
%                     maze.trials.types{end+1} = 'CL correct';
%                 elseif strcmp(GUI_states{this_evt(iT)}, 'C_L') && strcmp(GUI_states{this_evt(iT)+1}, 'R_Rew')
%                     maze.trials.types{end+1} = 'CL error';
%                     %                 end
%                 end
%             end
%         end
%         
%         if strcmp(states{ii}, 'F_R')
%             maze.trials.times = [maze.trials.times; maze.events.(states{ii})];
%             for iT = 1:length(this_evt)
%                 maze.trials.types{end+1} = 'FR correct';
%             end
%         end
%         
%         if strcmp(states{ii}, 'F_L')
%             maze.trials.times = [maze.trials.times; maze.events.(states{ii})];
%             for iT = 1:length(this_evt)
%                 maze.trials.types{end+1} = 'FL correct';
%             end
%         end
%         
%     else
%         this_evt = find(contains(GUI_states, states{ii}));
%         if this_evt(end) >= length(GUI_state_times)
%             this_evt = this_evt(1:end-1);
%         end
%         maze.events.(states{ii}) = GUI_state_times(this_evt+1);
%     end
% end
%% run the MAZE GUI to get User defined experimental time points.
% 
% 
% Maze_GUI_NLX_sandbox
% 
% 
% while ~exist('GUI_states', 'var')
% pause(0.1);
% end
% %% convert to trials

box_in = GUI_state_times(contains(GUI_states,  'Box_in'));
states = unique(GUI_states);

maze = [];
maze.trials.times = []; maze.trials.types = [];
for ii = length(states):-1:1
    
    if contains(states{ii}, {'Box_in', 'F', 'C'})
        this_evt = find(contains(GUI_states, states{ii}));
        if this_evt(end) == length(GUI_state_times)
            this_evt = this_evt(1:end-1);
        end
        
        maze.events.(states{ii}) = [GUI_state_times(this_evt);  GUI_state_times(this_evt+1)]';
        
        if strcmp(states{ii}, 'C_R') || strcmp(states{ii}, 'C_L')
            maze.trials.times = [maze.trials.times; maze.events.(states{ii})];
            for iT = 1:length(this_evt)
                if  strcmp(GUI_states{this_evt(iT)}, 'C_R') && strcmp(GUI_states{this_evt(iT)+1}, 'R_Rew')
                    maze.trials.types{end+1} = 'CR correct';
                elseif strcmp(GUI_states{this_evt(iT)}, 'C_R') && strcmp(GUI_states{this_evt(iT)+1}, 'L_Rew')
                    maze.trials.types{end+1} = 'CR error';
                elseif  strcmp(GUI_states{this_evt(iT)}, 'C_L') && strcmp(GUI_states{this_evt(iT)+1}, 'L_Rew')
                    maze.trials.types{end+1} = 'CL correct';
                elseif strcmp(GUI_states{this_evt(iT)}, 'C_L') && strcmp(GUI_states{this_evt(iT)+1}, 'R_Rew')
                    maze.trials.types{end+1} = 'CL error';
                    %                 end
                end
            end
        end
        
        if strcmp(states{ii}, 'F_R')
            maze.trials.times = [maze.trials.times; maze.events.(states{ii})];
            for iT = 1:length(this_evt)
                maze.trials.types{end+1} = 'FR correct';
            end
        end
        
        if strcmp(states{ii}, 'F_L')
            maze.trials.times = [maze.trials.times; maze.events.(states{ii})];
            for iT = 1:length(this_evt)
                maze.trials.types{end+1} = 'FL correct';
            end
        end
        
    else
        this_evt = find(contains(GUI_states, states{ii}));
        if this_evt(end) >= length(GUI_state_times)
            this_evt = this_evt(1:end-1);
        end
        maze.events.(states{ii}) = GUI_state_times(this_evt+1);
    end
end
maze.GUI.times = GUI_state_times;
maze.GUI.states = GUI_states;
% sort the trials by time

[~, idx] = sort(maze.trials.times(:,1));
maze.trials.times = maze.trials.times(idx,:);
maze.trials.types = maze.trials.types(idx);

save('maze.mat', 'maze', '-v7.3')
% 
%%  Check the events and split the trials into a summary plot.
figure(902)
clf
types = { 'FR correct'  'FL correct', 'CR correct', 'CR error', 'CL correct', 'CL error'};
c_ord = [0, 180, 120;...
    100 125 200;...
    0, 200, 100;...
    255, 0 , 0;...
    88 136 175;...
    255, 0 , 0]/255;
if length(maze.trials.times) >= 20
    n = 6;
    
else
    n= 4;
end
m = ceil(length(maze.trials.times)/n);

for ii  = 1:length(maze.trials.types)
    subplot(m,n,ii)
    hold on
    plot(pos.data(1,:), pos.data(2,:), '.', 'color', [.7 .7 .7 .2])
    
    plot(pos.data(1,nearest_idx3(maze.trials.times(ii,1), pos.tvec):nearest_idx3(maze.trials.times(ii,2), pos.tvec)), pos.data(2,nearest_idx3(maze.trials.times(ii,1), pos.tvec):nearest_idx3(maze.trials.times(ii,2), pos.tvec)), '.',...
        'color', c_ord(find(contains(types, maze.trials.types{ii})),:));
    ht = text(min(pos.data(1,:))*1.1, max(pos.data(2,:)*.98), [maze.trials.types{ii} '-' num2str(maze.trials.times(ii,1),6)]);
    ht = text(min(pos.data(1,:))*1.1, max(pos.data(2,:)*.90), ['         '  num2str(maze.trials.times(ii,2),6)]);
    
    xlim([min(pos.data(1,:)) max(pos.data(1,:))]);
    ylim([min(pos.data(2,:)) max(pos.data(2,:))]);
end
% 
% %% plot all events together
% % figure(903)
% % clf
% % c_ord = linspecer(length(maze.trials.types)+3);
% %
% % hold on
% % plot(pos.data(1,:), pos.data(2,:), '.', 'color', [.7 .7 .7 .2])
% %
% % for ii = size(maze.Box_in,1):-1:1
% %     plot(pos.data(1,nearest_idx3(maze.Box_in(ii,1), pos.tvec):nearest_idx3(maze.Box_in(ii,2), pos.tvec)), pos.data(2,nearest_idx3(maze.Box_in(ii,1), pos.tvec):nearest_idx3(maze.Box_in(ii,2), pos.tvec)), '.',...
% %         'color', c_ord(7,:));
% % end
% %
% % for ii = size(maze.F_L,1):-1:1
% %     plot(pos.data(1,nearest_idx3(maze.F_L(ii,1), pos.tvec):nearest_idx3(maze.F_L(ii,2), pos.tvec)), pos.data(2,nearest_idx3(maze.F_L(ii,1), pos.tvec):nearest_idx3(maze.F_L(ii,2), pos.tvec)), '.',...
% %         'color', c_ord(5,:));
% % end
% %
% % for ii = size(maze.F_R,1):-1:1
% %     plot(pos.data(1,nearest_idx3(maze.F_R(ii,1), pos.tvec):nearest_idx3(maze.F_R(ii,2), pos.tvec)), pos.data(2,nearest_idx3(maze.F_R(ii,1), pos.tvec):nearest_idx3(maze.F_R(ii,2), pos.tvec)), '.',...
% %         'color', c_ord(6,:));
% % end
% %
% % for ii  = 1:length(maze.trials.types)
% %         plot(pos.data(1,nearest_idx3(maze.trials.times(ii,1), pos.tvec):nearest_idx3(maze.trials.times(ii,2), pos.tvec)), pos.data(2,nearest_idx3(maze.trials.times(ii,1), pos.tvec):nearest_idx3(maze.trials.times(ii,2), pos.tvec)), '.',...
% %         'color', c_ord(find(contains(types, maze.trials.types{ii})),:));
% % end
% %
% % % legend(['Box', 'F_L', 'F_R', types])
% %
% %
% % xlim([min(pos.data(1,:)) max(pos.data(1,:))]);
% % ylim([min(pos.data(2,:)) max(pos.data(2,:))]);
% 
% 
% %% plot for DLC
% 
% %%  Check the events and split the trials into a summary plot.
% % figure(902)
% % clf
% % types = {'FR correct', 'FL correct', 'CR correct', 'CR error', 'CL correct', 'CL error'};
% % c_ord = [175, 255, 10;...
% %     170 75 170;...
% %     10, 186, 140;...
% %     200, 0 , 0;...
% %     75 125 200;...
% %     255, 0 , 0]/255;
% % 
% % for ii  = 1:length(maze.trials.types)
% %     subplot(4,4,ii)
% %     hold on
% %     plot(behav.position(:,1), behav.position(:,2), '.', 'color', [.7 .7 .7 .2])
% %     
% %     plot(behav.position(nearest_idx3(maze.trials.times(ii,1), behav.time):nearest_idx3(maze.trials.times(ii,2), behav.time),1), behav.position(nearest_idx3(maze.trials.times(ii,1), behav.time):nearest_idx3(maze.trials.times(ii,2), behav.time),2), '.',...
% %         'color', c_ord(find(contains(types, maze.trials.types{ii})),:));
% %     ht = text(min(behav.position(:,1))*1.1, max(behav.position(:,2)*.98), [maze.trials.types{ii} '-' num2str(maze.trials.times(ii,1),6)]);
% %     ht = text(min(behav.position(:,1))*1.1, max(behav.position(:,2)*.90), ['         '  num2str(maze.trials.times(ii,2),6)]);
% %     
% %     xlim([min(behav.position(:,1)) max(behav.position(:,1))]);
% %     ylim([min(behav.position(:,2)) max(behav.position(:,2))]);
% % end
