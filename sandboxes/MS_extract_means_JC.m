function [data_out, data_out_REM, data_out_SWS,Threshold,labels] = MS_extract_means_JC(cell_idx)
%% MS_extract_means_JC: used for pulling out the mean values for several measures and the time of the recording relative to the track experiment.
% The function will loop through all the session folders in the processed
% ms data and loads the ms_seg*.mat file.
%
% inputs
%     - cell_idx   [1 x n array] index for selecting only a subset of
%     cells. If empty it will use all the cells.
%
%    Outputs:
%    - data_out  [recording blocks x nVariables] will save these
%    variables, where iF is the recording block.
%        data_out(iF,1) = ms_seg.time2trk;
%        data_out(iF,2) = mean(sum(ms_seg.Binary/length(ms_seg.Binary)));  mean firing rate (Hz) across cells.
%        data_out(iF,3) = mean(ms_seg.d_amp);
%        data_out(iF,4) = mean(ms_seg.t_amp);
%        data_out(iF,5) = mean(ms_seg.LG_amp);
%        data_out(iF,6) = mean(ms_seg.HG_amp);
%        data_out(iF,7) = this_state;    % save the REM or SWS state. REM == 1 SWS ==0;
%
%
%
%
% EC 2020-08-19   initial version
%
%
%% load trk and update and extract time values.
if exist(['ms_trk.mat'], 'file')
    load('ms_trk.mat');
    
else
    error('no ms_trk.mat found')
end

% get time
time_str = regexp(ms_trk.file_names,'\d*','Match');
ms_trk.time_labels = datestr(datestr([time_str{1} ':' time_str{2} ':' time_str{3}]), 'HH:MM:SS');
temp_time = datetime(ms_trk.time_labels);
temp_time.Format = 'HH:mm:ss';
temp_time = temp_time+seconds((ms_trk.time(end)- ms_trk.time(1))/1000);

trk_time = datevec(ms_trk.time_labels);
trk_end_time = datevec(datestr(temp_time, 'HH:MM:SS'));

%% get the names of the folders in this dir and sort them based on date created.
this_dir = dir('*_M*');

% hack to sort based on recording time not file datenum.  
for iD = length(this_dir):-1:1
    fname_parts = strsplit(this_dir(iD).name, '_');
    temp_time = datetime([fname_parts{1}(2:end) ':' fname_parts{2}(2:end) ':' fname_parts{3}(2:end)]);
    temp_time.Format = 'HH:mm:ss'; 
    rec_time(iD) = temp_time; 
end
[~, rec_order_idx] = sort(rec_time, 'ascend'); 

this_dir = this_dir(rec_order_idx); 

% loop over each folder and get mean firing rate and mean amplitude
data_out_REM = [];
data_out_SWS = [];

labels = {'time2trk', 'binary', 'delta', 'theta', 'low gamma', 'mid gamma', 'high gamma', 'ultra high gamma', 'REM =1 SWS = 0'};
for iF = 1:length(this_dir)
    cd(this_dir(iF).name)
    f_load = FindFile_str(cd, 'resize');
    if isempty(f_load) || contains(this_dir(iF).name, 'remove') % if there is no ms_resize file skip this folder. 
        cd(this_dir(iF).folder)
        continue
    end
    load(f_load{1})
    if nargin == 0
        cell_idx = 1:size(ms_seg.Binary,2);
    end
    Fs = mode(diff(ms_seg.time));
    
    % check for time2trk
    if ~isfield(ms_seg, 'time2trk')
        seg_time = datevec(ms_seg.time_labels);
        % add time vector
        if strcmp(ms_seg.pre_post, 'pre')
            ms_seg.time2trk = etime(seg_time, trk_time)/60;
        elseif strcmp(ms_seg.pre_post, 'post')
            ms_seg.time2trk = etime(seg_time,trk_end_time)/60;
        end
    end

    if contains(this_dir(iF).name, 'REM')
        this_state = 1;
        data_out_REM(iF,1) = ms_seg.time2trk;
        data_out_REM(iF,2) = mean(sum(ms_seg.Binary(:,cell_idx))/(length(ms_seg.Binary)/Fs));
        data_out_REM(iF,3) = mean(ms_seg.d_amp);
        data_out_REM(iF,4) = mean(ms_seg.t_amp);
        data_out_REM(iF,5) = mean(ms_seg.LG_amp);
        data_out_REM(iF,6) = mean(ms_seg.MG_amp);
        data_out_REM(iF,7) = mean(ms_seg.HG_amp);
        data_out_REM(iF,8) = mean(ms_seg.UHG_amp);
    elseif contains(this_dir(iF).name, 'SWS') % maybe 'SW' just in case when i forgot to type SWS?
        this_state = 0;
        data_out_SWS(iF,1) = ms_seg.time2trk;
        data_out_SWS(iF,2) = mean(sum(ms_seg.Binary(:,cell_idx))/(length(ms_seg.Binary)/Fs));
        data_out_SWS(iF,3) = mean(ms_seg.d_amp);
        data_out_SWS(iF,4) = mean(ms_seg.t_amp);
        data_out_SWS(iF,5) = mean(ms_seg.LG_amp);
        data_out_SWS(iF,6) = mean(ms_seg.MG_amp);
        data_out_SWS(iF,7) = mean(ms_seg.HG_amp);
        data_out_SWS(iF,8) = mean(ms_seg.UHG_amp);
    end
    data_out(iF,1) = ms_seg.time2trk;
    data_out(iF,2) = mean(sum(ms_seg.Binary(:,cell_idx))/(length(ms_seg.Binary)/Fs));
    data_out(iF,3) = mean(ms_seg.d_amp);
    data_out(iF,4) = mean(ms_seg.t_amp);
    data_out(iF,5) = mean(ms_seg.LG_amp);
    data_out(iF,6) = mean(ms_seg.MG_amp);
    data_out(iF,7) = mean(ms_seg.HG_amp);
    data_out(iF,8) = mean(ms_seg.UHG_amp);
    % save the REM or SWS state. REM == 1 SWS ==0;
    data_out(iF,9) = this_state;
    
    %% Added by Jisoo
    
    Threshold = ms_seg.Binary_threshold;
    
    %% 
    cd(this_dir(iF).folder)
end

%% make a summary plot
REM_idx = find(data_out(:,9) == 1); % get REM indicies
SWS_idx = find(data_out(:,9)==0);
c_ord = linspecer(2); % set nice colours.

H=figure;
subplot(7,1,1)
hold on
plot(data_out(REM_idx,1), data_out(REM_idx, 2), '--*', 'color', c_ord(2,:))
plot(data_out(SWS_idx,1), data_out(SWS_idx, 2), '--*', 'color', c_ord(1,:))
legend('REM', 'SWS')
ylabel('Mean FR')

subplot(7,1,2)
hold on
plot(data_out(REM_idx,1), data_out(REM_idx, 3), '--*', 'color', c_ord(2,:))
plot(data_out(SWS_idx,1), data_out(SWS_idx, 3), '--*', 'color', c_ord(1,:))
legend('REM', 'SWS')
ylabel('Mean D')

subplot(7,1,3)
hold on
plot(data_out(REM_idx,1), data_out(REM_idx, 4), '--*', 'color', c_ord(2,:))
plot(data_out(SWS_idx,1), data_out(SWS_idx, 4), '--*', 'color', c_ord(1,:))
legend('REM', 'SWS')
ylabel('Mean T')

subplot(7,1,4)
hold on
plot(data_out(REM_idx,1), data_out(REM_idx, 5), '--*', 'color', c_ord(2,:))
plot(data_out(SWS_idx,1), data_out(SWS_idx, 5), '--*', 'color', c_ord(1,:))
legend('REM', 'SWS')
ylabel('Mean LG')

subplot(7,1,5)
hold on
plot(data_out(REM_idx,1), data_out(REM_idx, 6), '--*', 'color', c_ord(2,:))
plot(data_out(SWS_idx,1), data_out(SWS_idx, 6), '--*', 'color', c_ord(1,:))
legend('REM', 'SWS')
ylabel('Mean MG')

subplot(7,1,6)
hold on
plot(data_out(REM_idx,1), data_out(REM_idx, 7), '--*', 'color', c_ord(2,:))
plot(data_out(SWS_idx,1), data_out(SWS_idx, 7), '--*', 'color', c_ord(1,:))
legend('REM', 'SWS')
ylabel('Mean UG')

subplot(7,1,7)
hold on
plot(data_out(REM_idx,1), data_out(REM_idx, 8), '--*', 'color', c_ord(2,:))
plot(data_out(SWS_idx,1), data_out(SWS_idx, 8), '--*', 'color', c_ord(1,:))
legend('REM', 'SWS')
ylabel('Mean UHG')

saveas(H,'Dynamics_across_episodes.png') %added by jisoo

