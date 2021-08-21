function pREM_data = RB_get_pREM()
%% prep some data

c_ord = linspecer(5); % gen some nice colours. 1)blue 2)red 3)green 4)orange 4)purple

cfg_load = [];
lfp_chans = dir('*.ncs');

for ii = 1:length(lfp_chans)
    if strcmpi(lfp_chans(ii).name, 'CSC1.ncs')
       cfg_load.fc{ii} = []; 
    else
       cfg_load.fc{ii} = lfp_chans(ii).name; 
    end
end
cfg_load.fc(cellfun('isempty', cfg_load.fc)) = [];

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
legend({'infered hypnogram', 'scored data'})

%% load any spike files
if ~isempty(dir('*.t'))
    
    cfg = [];
    cfg.getTTnumbers = 0; % disable number loading.
    
    S  =  LoadSpikes(cfg);
else
    S = [];
end
%% get pREM

for iC = 1:length(csc.label)
   
    this_csc = csc;
    this_csc.data = csc.data(iC,:); 
    this_csc.label = csc.label{iC}; 
    this_csc.cfg.hdr = [];
    this_csc.cfg.hdr{1} = csc.cfg.hdr{iC}; 
    
    % congifuration
    cfg_pREM = [];
    cfg_pREM.min_len = 0.7;
    cfg_pREM.plot_flag = 1;
    cfg_pREM.REM_val = 3; 
    
    [pREM_idx{iC}, pREM_times{iC}, pREM_IV{iC}] = MS_get_pREM(this_csc, hypno == cfg_pREM.REM_val, cfg_pREM.min_len, [], cfg_pREM.plot_flag, S);
    
    h =  findobj('type','figure');
    
    for ih = 1:length(h)
        if h(ih).Number == 999 % fig number for the IPI plot. 
            saveas(h(ih), ['IPI_histogram_' this_csc.label '.png']);
        else
            saveas(h(ih), ['pREM_event_' num2str(h(ih).Number/1000) '_' this_csc.label(1:end-4) '.png']);
        end
    end
    close all

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

for iC = 1:length(csc.label)
    if isempty(pREM_times{iC})
        fprintf('<strong>0 pREM candidates detected on %s</strong>\n', csc.label{iC})
    else
       fprintf('<strong>%d pREM candidates detected on %s. Mean duration: %0.2f seconds</strong>\n', size(pREM_times{iC},1),csc.label{iC}, mean(pREM_times{iC}(:,2) - pREM_times{iC}(:,1)))
    end
end


% pick a channel and add in a raster to 
[~, best_chan] = max(cellfun('length', pREM_idx)); % use the csc with the most events. 


% collect the data for the output.  
parts = strsplit(cd, filesep);

pREM_data = [];
pREM_data.fname = parts{end}; 
pREM_data.csc = csc;
pREM_data.S = S;
pREM_data.hypno = hypno;
pREM_data.cfg = cfg_pREM; 
pREM_data.pREM.pREM_idx = pREM_idx;
pREM_data.pREM.pREM_times = pREM_times;
pREM_data.pREM.pREM_IV = pREM_IV;
pREM_data.best_chan = best_chan; 

save('pREM_data.mat', 'pREM_data', '-v7.3');
close all








