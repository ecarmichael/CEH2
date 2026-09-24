function states = HF_TFC_behav_cond(cfg_in, evts, type)
%% HF_TFC_Rec:  collects the occurance of events within a head fixed TFC 
% conditioning (typically day 3) for movement, behaviour, and any IV
% inputs.
%
% Inputs: cfg_in [struct]  contains configuration paramters. Any specified
%         will overwrite defaults. 
%
%         evts [struct]  OE events file.  Should contain all of the
%         movement, lick, tones, puffs, ....
%
%         iv_in [struct]  contains the start and stop times of some
%         interval data (typically SWR times]  (optional)

% plot the events as a check. 
figure;
%% initialize

cfg_def.start = '5'; 
cfg_def.tone1 = '11'; 
cfg_def.tone2 = '4'; 
cfg_def.puff = '13'; 

cfg_def.mov = '8'; 
cfg_def.lick = '14'; 

cfg = ProcessConfig(cfg_def, cfg_in);
%% get the behaviour states

switch type 


    case {'hab' 'Hab' 'TFC1' 'TFC2' 'Habituation' 'habituation'}
        states = [];
        states.task = iv(evts.t{ismember(evts.label, cfg.start)}');
        % context = [states.task.tstart(1) ;evts.t{ismember(evts.label, '4')}; states.task.tend(1)]; 
        states.contex_A = iv([0 870 1720]+ states.task.tstart , [360 1380 2060]+ states.task.tstart); 
        states.contex_B = iv([360 1380 2060]+ states.task.tstart, [[870 1720]+ states.task.tstart states.task.tend(1)]); 

        % states.contex_A = SelectIV([], states.contex, 1:2:length(states.contex.tstart)); 
        % states.contex_B = SelectIV([], states.contex, 2:2:length(states.contex.tstart)); 

        states.base_A = iv(evts.t{ismember(evts.label, '11')}(1:2:end)-20,evts.t{ismember(evts.label, '11')}(1:2:end));
        states.tone_A = iv(evts.t{ismember(evts.label, '11')}(1:2:end),evts.t{ismember(evts.label, '11')}(1:2:end)+20);
        states.trace_A = iv(evts.t{ismember(evts.label, '11')}(1:2:end)+20,evts.t{ismember(evts.label, '11')}(1:2:end)+35);

        states.base_B = iv(evts.t{ismember(evts.label, '13')}(1:2:end)-20,evts.t{ismember(evts.label, '13')}(1:2:end));
        states.tone_B = iv(evts.t{ismember(evts.label, '13')}(1:2:end),evts.t{ismember(evts.label, '13')}(1:2:end)+20);
        states.trace_B = iv(evts.t{ismember(evts.label, '13')}(1:2:end)+20,evts.t{ismember(evts.label, '13')}(1:2:end)+35);


        % split into Con A and B
        states.base_A_con_A = IntersectIV([], states.base_A, states.contex_A);
        states.tone_A_con_A = IntersectIV([], states.tone_A, states.contex_A);
        states.trace_A_con_A = IntersectIV([], states.trace_A, states.contex_A);

        states.base_A_con_B = IntersectIV([], states.base_A, states.contex_B);
        states.tone_A_con_B = IntersectIV([], states.tone_A, states.contex_B);
        states.trace_A_con_B = IntersectIV([], states.trace_A, states.contex_B);

        states.base_B_con_A = IntersectIV([], states.base_B, states.contex_A);
        states.tone_B_con_A = IntersectIV([], states.tone_B, states.contex_A);
        states.trace_B_con_A = IntersectIV([], states.trace_B, states.contex_A);

        states.base_B_con_B = IntersectIV([], states.base_B, states.contex_B);
        states.tone_B_con_B = IntersectIV([], states.tone_B, states.contex_B);
        states.trace_B_con_B = IntersectIV([], states.trace_B, states.contex_B);

% states = rmfield(states, {'contex'});

    case {'cond' 'Cond' 'TFC3', 'Conditioning'}
        states = [];
        states.task = iv(evts.t{ismember(evts.label, cfg.start)}');
        states.baseline = iv(evts.t{ismember(evts.label, cfg.tone1)}(1:2:end)-20,evts.t{ismember(evts.label, cfg.tone1)}(1:2:end));
        states.tone_A_con_A = iv(evts.t{ismember(evts.label, cfg.tone1)}(1:2:end),evts.t{ismember(evts.label, cfg.tone1)}(1:2:end)+20);
        states.trace_A_con_A = iv(evts.t{ismember(evts.label, cfg.tone1)}(1:2:end)+20, evts.t{ismember(evts.label, cfg.tone1)}(1:2:end)+35);
        states.US = iv(evts.t{ismember(evts.label, cfg.puff)}(1:2:end), evts.t{ismember(evts.label, cfg.puff)}(1:2:end)+2);


    case {'Rec' 'Test' 'TFC4' 'TFC5' 'Recall', 'rec' 'test' 'retrieval' 'Retrieval'}
        states = [];
        if isscalar(evts.t{ismember(evts.label, cfg.start)}')
            states.task = iv(evts.t{ismember(evts.label, cfg.start)}, evts.t{ismember(evts.label, '8')}(end));
        else
            states.task = iv(evts.t{ismember(evts.label, cfg.start)}');
        end
        states.base_A_con_B = iv(evts.t{ismember(evts.label, '11')}(1:2:end)-20,evts.t{ismember(evts.label, '11')}(1:2:end));
        states.tone_A_con_B = iv(evts.t{ismember(evts.label, '11')}(1:2:end),evts.t{ismember(evts.label, '11')}(1:2:end)+20);
        states.trace_A_con_B = iv(evts.t{ismember(evts.label, '11')}(1:2:end)+20, evts.t{ismember(evts.label, '11')}(1:2:end)+35);

        states.base_B_con_B = iv(evts.t{ismember(evts.label, '13')}(1:2:end)-20,evts.t{ismember(evts.label, '13')}(1:2:end));
        states.tone_B_con_B = iv(evts.t{ismember(evts.label, '13')}(1:2:end),evts.t{ismember(evts.label, '13')}(1:2:end)+20);
        states.trace_B_con_B = iv(evts.t{ismember(evts.label, '13')}(1:2:end)+20, evts.t{ismember(evts.label, '13')}(1:2:end)+35);

end



%% plot for carity

figure(5010)
clf

hold on
s_list = fieldnames(states);
c_ord = MS_linspecer(length(s_list));

for ii = 1:length(s_list)
    for jj = 1:length(states.(s_list{ii}).tstart)
        rectangle('position', [states.(s_list{ii}).tstart(jj), ii, (states.(s_list{ii}).tend(jj) - states.(s_list{ii}).tstart(jj)), 1], 'FaceColor', c_ord(ii,:)); 
    end
end

s_list = replace(s_list, '_', ' '); 
s_list = replace(s_list, '_', ' '); 

set(gca, 'YTick', 1.5:1:length(s_list)+.5, 'YTickLabel', s_list)