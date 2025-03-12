%% pola data loading and curration

% move to folder with ms.amt

% load it as a struct 

ms = load('ms.mat'); 

% rename deconv and denoise 

ms.denoise = ms.denoisedCa; 
ms.deconv = ms.deconvolvedSig; 

ms = rmfield(ms, 'denoisedCa'); 
ms = rmfield(ms, 'deconvolvedSig'); 


% add in the binarized trace and detrend data
ms = MS_append_filt_binary(ms);


% add in the detrend data from the filtered signal (makes for nicer plots)
ms.detrendRaw = []; % clear it just to be safe. 

    for trace_i = size(ms.FiltTraces,2):-1:1 % go backwards for speed. 
       ms.detrendRaw(:,trace_i)=detrend(ms.FiltTraces(:,trace_i),2);  
    end

% clean up the SFPs
ms = MS_append_sharp_SFPs(ms);

% get the time from the .csv file

ms = MS_append_timeStamps(ms, cd);

%% plot the data

MS_Ca_check(ms)


%% optional currate the data

% press 'space' to keep, 'delete' to remove', and 'back arrow' to redo a
% cell(s);

% simple version
[keep_idx, ms] = MS_curate_cells(ms); % appends the keep_idx to the ms struct. also gives it as an output; 

% version with ROI curration for speed. 
[keep_idx, ms] = MS_curate_cells(ms, [], 1); % appends the keep_idx to the ms struct. also gives it as an output; second input is for merging data. third is for using an ROI to remove cells outside the ROI. 

% remove the rejected cells

ms_out = MS_Ca_good_cells(ms); % we'll save another output struct just in case you need the full set. 


%% plot again

MS_Ca_check(ms_out)