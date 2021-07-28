%% pREM setup sandbox
% just a sandbox for testing the MS_get_pREM.m function. 


% RB data notes: Sample data for 1 mouse (220).
% Within the main folder are 2 experiment-specific subfolders (each mouse underwent testing under at least 2 different conditions (sometimes 3).
%AC is control data for mouse 220 (no MS GABAergic inhibition at any point of experiment), A is test group (MS GABAergic inhibition selectively during REMs).
% Within each experiment-specific subfolder is the nsc file from the best (largest, most stable recording) and the corresponding hypnogram. 
%(Each datapoint is 1 s, data was scored with 5 s resolution (non-overlapping windows); 
%1 = wake, 2 = NREMs, 3 = REMs, 4 = wake + MS GABAergic inhhibition, 5 = NREMs + MS GABAergic inhibition, 6 = REMs + MS GABAergic inhibition).


%% add codebases

addpath(genpath('C:\Users\ecarm\Documents\GitHub\CEH2'));
addpath(genpath('C:\Users\ecarm\Documents\GitHub\vandermeerlab\code-matlab\shared'));

% data_dir = 'J:\Scored_NSC_Data_files\Mouse 220\Experiment 1 (AC)';
% data_dir = 'J:\Scored_NSC_Data_files\Mouse 220\Experiment 2 (A)';
data_dir = 'J:\Williams_Lab\RB_data\220 AC'; 

cd(data_dir)
% LFP_dir = 'J:\Williams_Lab\Jisoo\LFP data\Jisoo'; 

c_ord = linspecer(5); % gen some nice colours. 1)blue 2)red 3)green 4)orange 4)purple

%% prep some data
cfg_load = [];
cfg_load.desired_sampling_frequency  = 1280; % closest to 1250 that FS of 32000 can get with whole decimation. 

csc = MS_LoadCSC(cfg_load);
Fs = csc.cfg.hdr{1}.SamplingFrequency; % get the sampling freq from nlx header. 
load('Hypnogram.mat')


%% fill in hypno to match csc length; hacky but works. FIX LATER
% hypno = round(interp1(1:length(Hypnogram), Hypnogram, csc.tvec));

hypno = NaN(size(csc.tvec));
for ii = 1:length(Hypnogram)
%     if (ii+Fs*ii) > length(csc.tvec)
%         break;
%     end
    if ii == 1
        hypno(ii:(Fs*ii)-1) = Hypnogram(ii);
        hypno(Fs*ii:(Fs*(ii+1))-1) = Hypnogram(ii+1);
    elseif ii == length(Hypnogram)
        hypno(Fs*ii:(Fs*(ii+1))-1) = Hypnogram(ii);
    else
        hypno(Fs*ii:(Fs*(ii+1))-1) = Hypnogram(ii+1);
    end
end

hypno = hypno(1:length(csc.tvec)); 
hypno(isnan(hypno)) = 5; % check for nans
% check 
figure(101)
hold on
this_tvec = csc.tvec-csc.tvec(1);
plot(this_tvec, hypno, 'color', c_ord(1,:)); % plot the new hypno. 
plot(0:length(Hypnogram)-1, Hypnogram, '--', 'color', c_ord(2,:)); % plot the original. 
xlim([0 length(Hypnogram)-1]); 
ylim([0 max(hypno)+1])

%% get pREM
REM_val = 3; % 

for iC = 1:length(csc.label)
    
    this_csc = csc;
    this_csc.data = csc.data(iC,:); 
    this_csc.label = csc.label{iC}; 
    this_csc.cfg.hdr = [];
    this_csc.cfg.hdr{1} = csc.cfg.hdr{iC}; 
    
    [pREM_idx{iC}, pREM_times{iC}, pREM_IV{iC}] = MS_get_pREM(this_csc, hypno == REM_val, 0.7, [], 1);
    
    h =  findobj('type','figure');
    
    for ih = 1:length(h)
        if ih == 1
            saveas(h(ih), ['ISI_histogram_' this_csc.label '.png']);
        else
            saveas(h(ih), ['pREM_event_' num2str(ih-1) '_' this_csc.label '.png']);
        end
    end
    close all

   fprintf('<strong>%d pREM candidates detected on %s. Mean duration: %0.2f seconds</strong>\n', length(pREM_times{iC}),csc.label{iC}, mean(pREM_times{iC}(:,2) - pREM_times{iC}(:,1)))
    %% add in the pREM times to the hypnogram plot.
%     figure(500+iC)
%     hold on
%     this_tvec = csc.tvec-csc.tvec(1);
%     plot(this_tvec, hypno, 'color', c_ord(1,:)); % plot the new hypno.
%     plot(0:length(Hypnogram)-1, Hypnogram, '--', 'color', c_ord(2,:)); % plot the original.
%     plot(this_tvec, (this_csc.data*100)+3.5,'k')
%     
%     xlim([0 length(Hypnogram)-1]);
%     ylim([0 max(hypno)+1])
%     for ii = 1:size(pREM_idx{iC},1)
%         xline(this_tvec(pREM_idx{iC}(ii,1)), '--k', 'Start');
%         xline(this_tvec(pREM_idx{iC}(ii,2)), '--m', 'End');
%     end
    
end
