function Evt = GetLightDurations(cfg_in, On, Off)



cfg_def =[];
cfg_def.round_fac = 4; % used to initially cope with timing that is not perfect.
cfg_def.min_nEvents = 6;  % minimum number of events per category.  
cfg = ProcessConfig2(cfg_def, cfg_in);

%% warnings
if length(On) ~= length(Off);
    error('Intervals are not of equal length')
end

%% get a histogram of all the light periods.  Use this to determine the different durations in the recording

d_light = round(Off-On,cfg.round_fac);
[N, light_dur] = hist(d_light, unique(d_light)); % rounding can be dangerous based on your light source timing
light_dur(N<cfg.min_nEvents) = [];

Evt.label = cell(1,length(light_dur));
Evt.ts = cell(1,length(light_dur));
for Ld = 1:length(light_dur)
    Evt.label{Ld} = num2str(light_dur(Ld));
    for iL = 1:length(On)
        if d_light(iL) == light_dur(Ld)
            Evt.ts{Ld} = cat(1,Evt.ts{Ld}, [On(iL), Off(iL)]);
        end
    end
end

%% put it all into the Evt format
Evt.type = 'ts';
Evt.cfg = cfg;

%% disp the outcome
fprintf(['\n' num2str(length(Evt.label)) ' Event(s) types found:\n']) 
    for iLabel = 1:length(Evt.label)
        n_spaces = 10 - length(Evt.label{iLabel});
        spacing_arg = ['%-', num2str(n_spaces),'s'];
        padded_string = sprintf(spacing_arg, ' ');
        fprintf([num2str(Evt.label{iLabel}) ':' padded_string   num2str(length(Evt.ts{iLabel})) '\n'])
    
    end
