%% check the timing

clearvars   -except behav_rec behav_enc Rad


%%

trl_s_idx = nearest_idx(Rad.D2025_02_20.JKA_07.recall.tstart, behav_rec.time);
trl_e_idx = nearest_idx(Rad.D2025_02_20.JKA_07.recall.tend, behav_rec.time);


display(1/mode(diff(behav_rec.time)))

figure(1919)
clf
hold on
plot(behav_rec.time, smooth(behav_rec.position(:,1),100))
plot(behav_rec.time, smooth(behav_rec.position(:,2),100))

hold on

for ii = 1:length(Rad.D2025_02_20.JKA_07.recall.tstart)
    xline(Rad.D2025_02_20.JKA_07.recall.tstart(ii), 'b', ['start: ' num2str(Rad.D2025_02_20.JKA_07.recall.tstart(ii))])% ' | '...
        %'idx time: ' num2str(behav_rec.time(trl_s_idx(ii)))])
    xline(Rad.D2025_02_20.JKA_07.recall.tend(ii), 'r', ['start: ' num2str(Rad.D2025_02_20.JKA_07.recall.tend(ii))]) %' | '...
        %'idx time: ' num2str(behav_rec.time(trl_e_idx(ii)))])
end
fprintf('Recording time: %2.1f | Log dur: %2.1f\n',behav_rec.time(end)- behav_rec.time(1),Rad.D2025_02_20.JKA_07.recall.tend(end))
