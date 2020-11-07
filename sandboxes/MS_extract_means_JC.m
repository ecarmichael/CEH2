function [data_out, data_out_REM, data_out_SWS,Threshold] = MS_extract_means_JC(cell_idx)
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
%
%% get the names of the folders in this dir and sort them based on date created.
this_dir = dir;
dirFlags = [this_dir.isdir];

this_dir = this_dir(dirFlags);
[~,idx] = sort([this_dir.datenum]);
this_dir = this_dir(idx);
keep_idx = zeros(1, length(this_dir));

for iD = 1:length(this_dir)
    if contains(this_dir(iD).name, '_M')
        keep_idx(iD) = 1;
    end
end

this_dir = this_dir(logical(keep_idx));

this_dir.name;

% loop over each folder and get mean firing rate and mean amplitude
data_out_REM = [];
data_out_SWS = [];
for iF = 1:length(this_dir)
    cd(this_dir(iF).name)
    f_load = FindFile_str(cd, 'ms_resize');
    if isempty(f_load) % if there is no ms_resize file skip this folder. 
        continue
    end
    load(f_load{1})
    if nargin == 0
        cell_idx = 1:size(ms_seg.Binary,2);
    end
    Fs = mode(diff(ms_seg.time));
    
    if contains(this_dir(iF).name, 'REM')
        this_state = 1;
        data_out_REM(iF,1) =  ms_seg.time2trk;
        data_out_REM(iF,2) = mean(sum(ms_seg.Binary(:,cell_idx))/(length(ms_seg.Binary)/Fs));
        data_out_REM(iF,3) = mean(ms_seg.d_amp);
        data_out_REM(iF,4) = mean(ms_seg.t_amp);
        data_out_REM(iF,5) = mean(ms_seg.LG_amp);
        data_out_REM(iF,6) = mean(ms_seg.HG_amp);
    elseif contains(this_dir(iF).name, 'SWS') % maybe 'SW' just in case when i forgot to type SWS?
        this_state = 0;
        data_out_SWS(iF,1) = ms_seg.time2trk;
        data_out_SWS(iF,2) = mean(sum(ms_seg.Binary(:,cell_idx))/(length(ms_seg.Binary)/Fs));
        data_out_SWS(iF,3) = mean(ms_seg.d_amp);
        data_out_SWS(iF,4) = mean(ms_seg.t_amp);
        data_out_SWS(iF,5) = mean(ms_seg.LG_amp);
        data_out_SWS(iF,6) = mean(ms_seg.HG_amp);
    end
    data_out(iF,1) = ms_seg.time2trk;
    data_out(iF,2) = mean(sum(ms_seg.Binary(:,cell_idx))/(length(ms_seg.Binary)/Fs));
    data_out(iF,3) = mean(ms_seg.d_amp);
    data_out(iF,4) = mean(ms_seg.t_amp);
    data_out(iF,5) = mean(ms_seg.LG_amp);
    data_out(iF,6) = mean(ms_seg.HG_amp);
    % save the REM or SWS state. REM == 1 SWS ==0;
    data_out(iF,7) = this_state;
    
    %% Added by Jisoo
    
    Threshold = ms_seg.Binary_threshold;
    
    %% 
    cd(this_dir(iF).folder)
end

%% make a summary plot
REM_idx = find(data_out(:,7) == 1); % get REM indicies
SWS_idx = find(data_out(:,7)==0);
c_ord = linspecer(2); % set nice colours.

H=figure;
subplot(5,1,1)
hold on
plot(data_out(REM_idx,1), data_out(REM_idx, 2), '--*', 'color', c_ord(2,:))
plot(data_out(SWS_idx,1), data_out(SWS_idx, 2), '--*', 'color', c_ord(1,:))
legend('REM', 'SWS')
ylabel('Mean FR')

subplot(5,1,2)
hold on
plot(data_out(REM_idx,1), data_out(REM_idx, 3), '--*', 'color', c_ord(2,:))
plot(data_out(SWS_idx,1), data_out(SWS_idx, 3), '--*', 'color', c_ord(1,:))
legend('REM', 'SWS')
ylabel('Mean D')

subplot(5,1,3)
hold on
plot(data_out(REM_idx,1), data_out(REM_idx, 4), '--*', 'color', c_ord(2,:))
plot(data_out(SWS_idx,1), data_out(SWS_idx, 4), '--*', 'color', c_ord(1,:))
legend('REM', 'SWS')
ylabel('Mean T')

subplot(5,1,4)
hold on
plot(data_out(REM_idx,1), data_out(REM_idx, 5), '--*', 'color', c_ord(2,:))
plot(data_out(SWS_idx,1), data_out(SWS_idx, 5), '--*', 'color', c_ord(1,:))
legend('REM', 'SWS')
ylabel('Mean LG')

subplot(5,1,5)
hold on
plot(data_out(REM_idx,1), data_out(REM_idx, 6), '--*', 'color', c_ord(2,:))
plot(data_out(SWS_idx,1), data_out(SWS_idx, 6), '--*', 'color', c_ord(1,:))
legend('REM', 'SWS')
ylabel('Mean HG')
saveas(H,'Dynamics_across_episodes.png') %added by jisoo

