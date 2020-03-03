function data_ft = MS_TSDtoFT(cfg_in,data)
% function data_ft = TSDtoFT(cfg_in,data_tsd)
%
% converts tsd into fieldtrip (ft) data structure
%
% cfg_def.mode = 'as-is'; % {'as-is','resample'}, defines how to deal with gaps in data
%
% MvdM 2014-11-12 initial version
% NOTE multiple channels not yet implemented!
% EC modified version to avoid conflicts.
%       - fills in missing hdr information

cfg_def = [];
cfg_def.mode = 'as-is'; % {'as-is','resample'}, defines how to deal with gaps in data
cfg = ProcessConfig(cfg_def,cfg_in);

if ~CheckTSD(data)
    return;
end

dts = unique(diff(data.tvec));
if length(dts) > 1
    fprintf('\nTSDtoFT.m: WARNING: tvec diffs are not constant, cannot determine Fs.');
    % could approximate with a median if matches average closely enough
    Fs = 1./median(dts);
    fprintf('\nTSDtoFT.m: Fs %.2f estimated.\n',Fs);
else
    Fs = 1./dts;
    fprintf('\nTSDtoFT.m: Fs %.2f detected.\n',Fs);
end

%
data_ft = [];
switch cfg.mode
    
    case 'as-is'
        nSamples = size(data.data,2);
        
        data_ft.trial{1}  = data.data;
        data_ft.time{1}   = data.tvec;
        
    case 'resample' % currently does interpolation of data, could be improved with options like inserting NaNs or zeros if no sample nearby
        
        data_ft.time{1} = data.tvec(1):1./Fs:data.tvec(end);
        data_ft.trial{1} = interp1(data.tvec,data.data,data_ft.time{1},'nearest');
        
        nSamples = length(data_ft.time{1});
    otherwise
        error('Unknown mode.');
end

for ii = 1:size(data.data,1)
    ts    = data.tvec;
    mn(ii) = ts(1);
    mx(ii) = ts(end);
    ts1   = ts(1);
    md(ii) = mode(diff(double(ts-ts1)));
    
    % get the minimum and maximum across all channels
    if ii>1 && mn(ii)<min_all
        min_all = mn(ii);
    else
        min_all = mn(ii);
    end
    
    if ii>1 && mx(ii)>max_all
        max_all = mx(ii);
    else
        max_all = mx(ii);
    end
end

mode_dts  = mode(md);
rng       = double(max_all-min_all); % this is small num, can be double
tsinterp  = [0:mode_dts:rng]; % the timestamp interpolation axis

if size(data_ft.trial{1},1) ~= 1
    data_ft.trial{1} = data_ft.trial{1}';
end

if size(data_ft.time{1},1) ~= 1
    data_ft.time{1} = data_ft.time{1}';
end

% define hdr to match ft format.
data_ft.hdr.Fs = Fs;
data_ft.hdr.label = data.label; 
data_ft.hdr.nChan = size(data.data,1);
data_ft.hdr.nTrials = size(data.data,1); 
data_ft.hdr.nSamplesPre = 0; 
data_ft.hdr.nSamples = nSamples;
data_ft.hdr.TimeStampPerSample = mode_dts;
data_ft.hdr.FirstTimeStamp     = min_all;
data_ft.hdr.chantype = {'unknown'};
data_ft.hdr.chanunit = {data.units}; 
data_ft.hdr.LastTimeStamp      = uint64(tsinterp(end)) + min_all;


data_ft.label   = data.label;
data_ft.sampleinfo = [1 data_ft.hdr.nSamples];
data_ft.fsample = Fs;
data_ft.cfg = cfg;
end