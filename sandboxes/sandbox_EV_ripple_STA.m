
all_ripple = [];
for iCell = 1:length(ripple_lfp)
    
   all_ripple = [all_ripple, ripple_lfp{iCell}]; 
end

all_ripple = ripple_lfp{1};
%% 
tvec = (0:length(all_ripple))/Target_Sampling_Frequency;
tvec = tvec(1:end-1);
sta = mean(all_ripple,2);
sta_sem = std(all_ripple, [], 2)/sqrt(length(all_ripple)); 

figure
hold on
plot(tvec, sta, 'k'); 
plot(tvec, sta+sta_sem, '--r'); 
plot(tvec, sta-sta_sem, 'k'); 


%% spectrogram

% Triggered_Spec_FT(