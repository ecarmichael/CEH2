function trials = MS_NLX_append_MAZE(pos, evt_r)
%%  MS_NLX_append_MAZE:
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
% EC 2022-08-03   initial version based on MS_behav_append_MAZE
%
%
%
%% initialize

    
if  ~exist('pos', 'var')
    if exist('VT1.nvt', 'file')
        pos = LoadPos([]); 
    else
        error('No pos input or VT1.nvt found')
    end
end
    
if length(evt_r.t{find(contains(evt_r.label, 'Starting Recording'))}) > 1
    error('Events not restricted to time on maze')
end
  

pos_out = pos; 

%% write the event times for each trial and the type of trial.




% get all the trial starts
trials.choice = evt_r.t{strcmpi(evt_r.label, 'c')};
trials.start_C = evt_r.t{contains(evt_r.label, 'Box Exited')}((nearest_idx3(trials.choice, evt_r.t{contains(evt_r.label, 'Box Exited')},-1)));

trials.force_L = evt_r.t{strcmpi(evt_r.label, 'fl')};
trials.start_FL = evt_r.t{contains(evt_r.label, 'Box Exited')}(nearest_idx3(trials.force_L, evt_r.t{contains(evt_r.label, 'Box Exited')},-1));

trials.force_R = evt_r.t{strcmpi(evt_r.label, 'fr')};
trials.start_FR = evt_r.t{contains(evt_r.label, 'Box Exited')}((nearest_idx3(trials.force_R, evt_r.t{contains(evt_r.label, 'Box Exited')},-1)));

trials.rew_L = evt_r.t{strcmpi(evt_r.label, 'L')};
trials.enter_L = evt_r.t{contains(evt_r.label, 'L reward Entered')}((nearest_idx3(trials.rew_L, evt_r.t{contains(evt_r.label, 'L reward Entered')})));
% catch cases where the zoned video failed and default to UI value
if size(trials.enter_L, 1) > size(trials.enter_L,2)
    trials.enter_L = trials.enter_L'; 
end
L_flags = find(abs(trials.rew_L - trials.enter_L) > 30);
trials.enter_L(L_flags) = trials.rew_L(L_flags); 


trials.rew_R = evt_r.t{strcmpi(evt_r.label, 'R')};
trials.enter_R = evt_r.t{contains(evt_r.label, 'R reward Entered')}((nearest_idx3(trials.rew_R, evt_r.t{contains(evt_r.label, 'R reward Entered')})));
% catch cases where the zoned video failed and default to UI value
% catch cases where the zoned video failed and default to UI value
if size(trials.enter_R, 1) > size(trials.enter_R,2)
    trials.enter_R = trials.enter_R'; 
end
R_flags = find(abs(trials.rew_R - trials.enter_R) > 30);
trials.enter_R(R_flags) = trials.rew_R(R_flags); 

%     trials.box = evt_r.t{strcmpi(evt_r.label, 'b')};
%     trials.enter_box = evt_r.t{contains(evt_r.label, 'Box Entered')}((nearest_idx3(trials.box , evt_r.t{contains(evt_r.label, 'Box Entered')}, -1)));

trials.box = sort([trials.start_C, trials.start_FL,trials.start_FR]) - 30;

% split the choice trials into L / R 
type_vec = [repmat(1, 1, length(trials.start_C)), repmat(3, 1, length(trials.start_FR)), repmat(4, 1, length(trials.start_FL))]; 
time_vec = [trials.start_C, trials.start_FR, trials.start_FL]; 
[time_vec_s, idx] = sort(time_vec); 
type_vec = type_vec(idx); 

% types: C = 1, FR = 3, FL = 4, CR = 13, CL = 14; CR-C = 131, CR-I = 130,
% CL-C = 141, CL-I = 140; 
trials.start_CR = nan(1,length(type_vec)); trials.start_CL = trials.start_CR; 
for ii = 1:length(type_vec)
    if type_vec(ii) == 3 || type_vec(ii) == 4
        continue
    else
        if type_vec(ii-1) == 3
            trials.start_CL(ii) = time_vec_s(ii);
            fprintf('CL at %d, previous %d\n', time_vec_s(ii), type_vec(ii-1))
            type_vec(ii) = 14;
        elseif type_vec(ii-1) == 4
            trials.start_CR(ii) = time_vec_s(ii);
            fprintf('CR at %d, previous %d\n', time_vec_s(ii), type_vec(ii-1))
            type_vec(ii) = 13;

        end
    end
end
trials.start_CL(isnan(trials.start_CL)) = []; % remove empty spots. 
trials.start_CR(isnan(trials.start_CR)) = []; 


% determine if the trial was correct or incorrect based on previous
% outcome. 
% Choice L correct
type_vec_cl = []; type_vec_cr = [];
correct_cl = []; correct_cr = []; 
for ii = 1:length(trials.start_CL)
    R_diff = trials.enter_R - trials.start_CL(ii);
    L_diff = trials.enter_L - trials.start_CL(ii);
    fprintf('Trial %.0f: time to R %.0f | time to L %.0f  = ',ii, min(R_diff(R_diff> 0)), min(L_diff(L_diff> 0)))
    if (isempty(L_diff(L_diff> 0)) ||  min(R_diff(R_diff> 0)) < min(L_diff(L_diff> 0)))
        fprintf('R rew -> CL - C\n')
        type_vec_cl(ii) = 141;
        correct_cl(ii) = 1;
    elseif (isempty(R_diff(R_diff> 0)) ||  min(R_diff(R_diff> 0)) > min(L_diff(L_diff> 0)))
        fprintf('L rew -> CL - I\n')
        type_vec_cl(ii) = 140;
        correct_cl(ii) = 0;
    end
end

%Choice R correct
for ii = 1:length(trials.start_CR)
    R_diff = trials.enter_R - trials.start_CR(ii);
    L_diff = trials.enter_L - trials.start_CR(ii);
    fprintf('Trial %.0f: time to R %.0f | time to L %.0f  = ',ii, min(R_diff(R_diff> 0)), min(L_diff(L_diff> 0)))
    if (isempty(L_diff(L_diff> 0)) ||  min(R_diff(R_diff> 0)) < min(L_diff(L_diff> 0)))
        fprintf('R rew -> CR - I\n')
        type_vec_cr(ii) = 130;
        correct_cr(ii) = 0;
    elseif (isempty(R_diff(R_diff> 0)) ||  min(R_diff(R_diff> 0)) > min(L_diff(L_diff> 0)))
        fprintf('L rew -> CR - C\n')
        type_vec_cr(ii) = 131;
        correct_cr(ii) = 1;
    end
end


type_vec(find(type_vec == 14)) = type_vec_cl; 
type_vec(find(type_vec == 13)) = type_vec_cr; 
correct_trials = nan(size(type_vec));
correct_trials(find(type_vec == 141 | type_vec == 140)) = correct_cl; 
correct_trials(find(type_vec == 131 | type_vec == 130)) = correct_cr; 

type = {};
for ii = length(type_vec):-1:1
    if type_vec(ii) ==3
        type{ii} = 'FL';
    elseif type_vec(ii) ==4
        type{ii} = 'FR';
    elseif type_vec(ii) ==141
        type{ii} = 'CL-C';
    elseif type_vec(ii) ==140
        type{ii} = 'CL-I';
    elseif type_vec(ii) ==131
        type{ii} = 'CR-C';
    elseif type_vec(ii) ==130
        type{ii} = 'CR-I';
    end
end


% trials.start_CR = []; trials.start_CL = []; 
% for ii = 1:length(trials.start_C)
%     R_diff = trials.enter_R - trials.start_C(ii);
%     L_diff = trials.enter_L - trials.start_C(ii);
%     fprintf('Trial %.0f: time to R %.0f | time to L %.0f  = ',ii, min(R_diff(R_diff> 0)), min(L_diff(L_diff> 0)))
%     
%     % see which reward is the closest to the choice point. 
%     if  ~isempty(R_diff(R_diff> 0)) && (isempty(L_diff(L_diff> 0)) ||  min(R_diff(R_diff> 0)) < min(L_diff(L_diff> 0)))
%         trials.start_CR = [ trials.start_CR trials.start_C(ii)];
%         fprintf('R trial\n')
%     elseif ~isempty(L_diff(L_diff> 0)) && (isempty(R_diff(R_diff> 0)) || min(R_diff(R_diff> 0)) > min(L_diff(L_diff> 0)))
%         trials.start_CL = [ trials.start_CL trials.start_C(ii)];
%         fprintf('L trial\n')
%     end
% end


% split into correct incorrect 

% %vectorize trial times/types
% type_vec = [repmat(1, 1, length(trials.start_CL)), repmat(2, 1, length(trials.start_CR)), repmat(3, 1, length(trials.start_FL)), repmat(4, 1, length(trials.start_FR))]; 
% time_vec = [trials.start_CL, trials.start_CR, trials.start_FL, trials.start_FR]; 
% 
% [time_vec_s, idx] = sort(time_vec); 
% type_vec = type_vec(idx); 
% 
% correct = NaN(size(trials.start_C)); 
% type = {}; 
% for ii = length(type_vec):-1:1
%     if type_vec(ii) == 3
%         correct(ii) = NaN;
%         type{ii} = 'FL';
%         continue
%     elseif type_vec(ii) == 4
%         correct(ii) = NaN;
%         type{ii} = 'FR';
%         continue
%     end
%     
%     if type_vec(ii) == 1 && type_vec(ii-1) == 4
%         correct(ii) = 1; 
%         type{ii} = 'CL - C'; 
%     elseif type_vec(ii) == 1 && type_vec(ii-1) == 3
%         correct(ii) = 0; 
%         type{ii} = 'CL - I'; 
%     elseif type_vec(ii) == 2 && type_vec(ii-1) == 4
%         correct(ii) = 0; 
%         type{ii} = 'CR - I'; 
%     elseif type_vec(ii) == 2 && type_vec(ii-1) == 3
%         correct(ii) = 1; 
%         type{ii} = 'CR - C'; 
%     end
% end


trials.tstart = sort([trials.start_C, trials.start_FL,trials.start_FR]);
trials.tend = sort([trials.rew_L, trials.rew_R]) +3;

if length(trials.tstart) ~= length(trials.tend)
    if length(trials.tstart)-1 == length(trials.tend)
        warning('tstart (%.0f) tend (%.0f) differ by 1. Assuming incomplete last trial', length(trials.tstart),length(trials.tend));
    else
        error('%s: tstart (%.0f) tend (%.0f) differ. Determine why...',mfilename, length(trials.tstart),length(trials.tend));
        
    end
end



%% add in the correct/incorrect to both
trials.correct = correct_trials;
% same for the types
trials.type = type;
% collect outputs
pos_out.trials_nlx = trials;

%% test event
%
% % close(101)
% figure(101)
% subplot(1,2,1)
% pos_r = restrict(pos, evt.t{1}(this_rec), evt.t{2}(this_rec));
% hold on
% plot(pos_r.data(1,:), pos_r.data(2,:), '.k')
% 
% plot(pos_r.data(1,nearest_idx3(trials.start_C, pos_r.tvec)), pos_r.data(2,nearest_idx3(trials.start_C, pos_r.tvec)), 'sb', 'markersize', 12, 'markerfacecolor', 'b')
% plot(pos_r.data(1,nearest_idx3(trials.start_FR, pos_r.tvec)), pos_r.data(2,nearest_idx3(trials.start_FR, pos_r.tvec)), 'sr', 'markersize', 12, 'markerfacecolor', 'r')
% plot(pos_r.data(1,nearest_idx3(trials.start_FL, pos_r.tvec)), pos_r.data(2,nearest_idx3(trials.start_FL, pos_r.tvec)), 'sm', 'markersize', 12, 'markerfacecolor', 'm')
% 
% plot(pos_r.data(1,nearest_idx3(trials.enter_R, pos_r.tvec)), pos_r.data(2,nearest_idx3(trials.enter_R, pos_r.tvec)), 'dr', 'markersize', 12, 'markerfacecolor', 'r')
% plot(pos_r.data(1,nearest_idx3(trials.enter_L, pos_r.tvec)), pos_r.data(2,nearest_idx3(trials.enter_L, pos_r.tvec)), 'db', 'markersize', 12, 'markerfacecolor', 'b')
% 
% 
% plot(pos_r.data(1,nearest_idx3(trials.box, pos_r.tvec)), pos_r.data(2,nearest_idx3(trials.box, pos_r.tvec)), 'ob', 'markersize', 12, 'markerfacecolor', 'b')
% 
% % plot(pos_r.data(1,nearest_idx3(evt.t{find(contains(evt.label, 'Box Exited'))}, pos_r.tvec)), pos_r.data(2,nearest_idx3(evt.t{find(contains(evt.label, 'Box Exited'))}, pos_r.tvec)), 'sb', 'markersize', 12, 'markerfacecolor', 'b')
% % plot(pos_r.data(1,nearest_idx3(evt.t{find(contains(evt.label, 'fl'))}, pos_r.tvec)), pos_r.data(2,nearest_idx3(evt.t{find(contains(evt.label, 'fl'))}, pos_r.tvec)), 'sr', 'markersize', 12, 'markerfacecolor', 'r')
% 
% legend({'tracking', 'choise', 'Force R', 'Force L', 'Rew R', 'Rew L', 'Box ent'})
% 
% subplot(1,2,2)
% c_ord = linspecer(length(trials.tend));
% hold on
% for ii  = 1:length(trials.tend)
%     this_pos = restrict(pos, trials.tstart(ii), trials.tend(ii));
% 
%     plot(this_pos.data(1,:), this_pos.data(2,:), '.', 'color', c_ord(ii,:));
%     pause(0.5)
% end
%

%%  event plots using behav.
figure(903)
% clf
% subplot(1,2,2); 

c_ord = linspecer(length(trials.tend));

for ii  = 1:length(trials.tend)
    subplot(4,4,ii)
    hold on
plot(pos.data(1,:), pos.data(2,:), '.', 'color', [.7 .7 .7 .2])

%     this_pos = nearest_idx3(trials_ms.tstart(ii), behav.time);
    plot(pos.data(1,nearest_idx3(trials.tstart(ii), pos.tvec):nearest_idx3(trials.tend(ii), pos.tvec)), pos.data(2,nearest_idx3(trials.tstart(ii), pos.tvec):nearest_idx3(trials.tend(ii), pos.tvec)), '.', 'color', c_ord(ii,:));
    ht = text(min(pos.data(1,:))*1.1, max(pos.data(2,:)*.95), [trials.type{ii} '-' num2str(trials.tstart(ii),6)]);
    xlim([min(pos.data(1,:)) max(pos.data(1,:))]); 
ylim([min(pos.data(2,:)) max(pos.data(2,:))]);
%     drawnow
%         delete(ht); 

end
%     pause(1)

% close(903)
%%  tracking playback (for debugging)
figure(10101)
clf

plot(pos.data(1,:), pos.data(2,:), '.', 'color', [.7 .7 .7 .2])
hold on

for ii = 1:5:length(pos.tvec)
   hxy = plot(pos.data(1,ii), pos.data(2,ii), 'or'); 
   ht = text(min(pos.data(1,:))*1.1, max(pos.data(2,:)*.9), num2str(pos.tvec(ii)));
   drawnow
   pause(0.01)
   delete(hxy); 
   delete(ht)
end

