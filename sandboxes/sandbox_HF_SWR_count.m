% sandbox_HF_SWR_count




%% load some data

cfg_csc = []; 
% cfg_csc.desired_sampling_frequency = 2000; 
cfg_csc.fc = {'CSC52.ncs'}; 


csc = LoadCSC(cfg_csc);


evts = MS_LoadEvents

%% plot

figure(1)
clf

hold on
plot(csc.tvec - csc.tvec(1), csc.data,'k');

for ii = 1:length(evts.t{3})
    if ii == 10
xline(evts.t{3}(ii) - csc.tvec(1), '-b')
    else
    xline(evts.t{3}(ii) - csc.tvec(1), '-r')
    
    end

end
for ii = 1:length(evts.t{1})

    xline(evts.t{1}(ii) - csc.tvec(1), '--m')
end

%% restrict to long on off blue

csc_swr = restrict(csc, evts.t{3}(1),evts.t{4}(1)) 

%% detect SWR


MS_SWR_detector(csc_swr, csc.label{1}, 1)


