function csc_tsd = MS_LoadCSC(cfg_in)
% function csc_tsd = LoadCSC(cfg)
%
% loads Neuralynx .ncs files, handling missing samples correctly
%
% INPUTS: only cfg
%
% cfg.fc: cell array containing filenames to load
%   if no file_list field is specified, loads all *.Ncs files in current dir
% cfg.TimeConvFactor = 10^-6; % from nlx units to seconds
% cfg.VoltageConvFactor = 1; % if 0, return raw AD bits. If 1, data returned is in V
% cfg.decimateByFactor = []; % if non-empty, decimate data by specified factor
% cfg_def.desired_sampling_frequency = [];  % put in what you want the final FS to be.  Commonly 2kHz
% cfg.verbose = 1; Allow or suppress displaying of command window text
%
% OUTPUTS:
%
% csc_tsd: CSC data tsd struct
%
% MvdM 2014-06-18, 25 (use cfg_in)
% youkitan 2016-02-20 (fixed reading Header file)
% EC 2020-01-06 (made this more universal.  Will make own version in
% future.)
%

cfg_def.fc = {};
cfg_def.label = {};
cfg_def.TimeConvFactor = 10^-6; % 10^-6 means convert nlx units to seconds
cfg_def.VoltageConvFactor = 1; % 1 means output in volts, 1000 in mV, 10^6 in uV
cfg_def.decimateByFactor = []; % how much would you like to
cfg_def.desired_sampling_frequency = [];  % put in what you want the final FS to be.  Commonly 2kHz
cfg_def.verbose = 1;

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun); % this takes fields from cfg_in and puts them into cfg

% check for an attempt to both automate and mannually decimate
if ~isempty(cfg.decimateByFactor) && ~isempty(cfg.desired_sampling_frequency)
    error('Looks like you are trying to decimate the data both mannually (cfg.decimateByFactor) and automatically (cfg.desired_sampling_frequency). Use one or the other. Leave the undesied method empty')
end

if isempty(cfg.fc) % no filelist provided, load everything
    
    cfg.fc = FindFiles('*.ncs'); % swapped to lower case.
    
else
    
    if ~isa(cfg.fc,'cell')
        error('LoadCSC: cfg.fc should be a cell array.');
    end
    
end

if isempty(FindFiles(cfg.fc{1}))
    error('File does not exist in directory. Check spelling.');
end

fc = sort(cfg.fc);
nFiles = length(fc);

% track sample counts for different files, throw error if not equal
sample_count_tvec = nan(nFiles,1);
sample_count_data = nan(nFiles,1);

csc_tsd = tsd;

if cfg.VoltageConvFactor == 1
    csc_tsd.units = 'V';
elseif cfg.VoltageConvFactor == 1000
    csc_tsd.units = 'mV';
elseif cfg.VoltageConvFactor == 10^6
    csc_tsd.units = 'uV';
elseif cfg.VoltageConvFactor == 0
    csc_tsd.units = 'ADbits';
else
    error('Input voltage conversion factor is not a valid value.')
end

if cfg.verbose; fprintf('%s: Loading %d file(s)...\n',mfun,nFiles); end

for iF = 1:nFiles
    
    fname = fc{iF};
    
    % load raw data
    [Timestamps, ~, SampleFrequencies, NumberOfValidSamples, Samples, Header] = Nlx2MatCSC(fname, [1 1 1 1 1], 1, 1, []);
    
    % disabled channels cannot be loaded
    if Timestamps == 0
        error(['No csc data (disabled tetrode channel). Consider deleting ',fname,'.']);
    end
    
    % check for each data
    
    % check for constant sampling frequency
    Fs = unique(SampleFrequencies);
    if length(Fs) ~= 1
        error('More than one sampling frequency found!');
    end
    
    % extract information from header
    hdr = readCSCHeader(Header);
    
    % data format is n columns (blocks) of 512 samples, with a timestamp for each first
    % sample per block
    %
    % b1_1    b2_1    b3_1    ...   bn_1
    % b1_2    b2_2    b3_2    ...   bn_2
    % .       .       .       ...   .
    % .       .       .       ...   .
    % b1_512  b2_512  b3_512  ...   bn_512
    %
    % but, not all blocks have 512 good samples, some have less
    % (number is given in NumberOfValidSamples)
    
    % convert units
    Timestamps = Timestamps .* cfg.TimeConvFactor;
    if cfg.VoltageConvFactor ~= 0
        Samples = Samples .* cfg.VoltageConvFactor .* hdr.ADBitVolts;
    end
    
    % construct within-block tvec
    nSamplesPerBlock = size(Samples,1);
    block_tvec = 0:1./Fs:(nSamplesPerBlock-1)./Fs;
    
    % slow but tractable: loop over blocks remembering to initialize variables
    % first
    
    data = nan(numel(Samples),1); % allocate memory and then fill; can trim later
    tvec = nan(numel(Samples),1);
    
    idx = 1; % move this along as we go through the loop over all samples
    
    nBlocks = length(Timestamps);
    badBlocks = 0; % counter
    for iB = 1:nBlocks
        
        nvs = NumberOfValidSamples(iB);
        if nvs ~= 512, badBlocks = badBlocks + 1; end
        
        currentData = Samples(1:nvs,iB);
        currentTime = block_tvec(1:nvs)+Timestamps(iB);
        
        data(idx:idx+nvs-1) = currentData;
        tvec(idx:idx+nvs-1) = currentTime;
        
        idx = idx + nvs;
        
    end % of block loop
    
    cfg.badBlocks = badBlocks;
    
    if cfg.verbose; fprintf('%s: %s %d/%d bad blocks found (%.2f%%).\n',mfun,fname,badBlocks,nBlocks,(badBlocks./nBlocks).*100); end
    
    
    % remove nans
    data = data(~isnan(data));
    tvec = tvec(~isnan(tvec));
    
    % track sizes
    sample_count_tvec(iF) = length(tvec);
    sample_count_data(iF) = length(data);
    
    
    % decimate data if specified
    if ~isempty(cfg.decimateByFactor) && isempty(cfg.desired_sampling_frequency) % make sure you are not trying to decimate and not also automating the downsampling.
        
        fprintf('%s: Decimating by factor %d...\n',mfun,cfg.decimateByFactor)
        data = decimate(data,cfg.decimateByFactor);
        tvec = tvec(1:cfg.decimateByFactor:end);
        hdr.SamplingFrequency = hdr.SamplingFrequency./cfg.decimateByFactor;
        
    end
    
    % "automated decimation" (get a desired sampling frequency given the
    % data). Save the time/memory from having to do this after.
    % decimate data if specified
    if ~isempty(cfg.desired_sampling_frequency) && isempty(cfg.decimateByFactor) % make sure you are not trying to decimate and not also automating the downsampling.
        
        this_Fs = hdr.SamplingFrequency;
        this_decimate = round(this_Fs / cfg.desired_sampling_frequency); % must be a whole number.
        fprintf('%s: Current Sampling Frequency: %dHz...\n',mfun,this_Fs);
        fprintf('%s: Decimating by factor %d to output data at desired Fs of %dHz \n',mfun,this_decimate, cfg.desired_sampling_frequency)
        
        
        data = decimate(data,this_decimate);
        tvec = tvec(1:this_decimate:end);
        hdr.SamplingFrequency = hdr.SamplingFrequency./this_decimate;
        
    end
    
        % check if the data is the same length for each channel. EC moved
        % this so that it can actually decimate. 
    if iF >1 && length(data) ~= length(csc_tsd.data(iF-1,:))
        message = 'Data lengths differ across channels.';
        error(message);
    end
    
    
    % done, add to tsd
    csc_tsd.tvec = tvec;
    csc_tsd.data(iF,:) = data;
    
    if isempty(cfg.label)
        [~,fn,fe] = fileparts(fname);
        csc_tsd.label{iF} = cat(2,fn,fe);
    else
        csc_tsd.label{iF} = cfg.label{iF};
    end
    
    csc_tsd.cfg.hdr{iF} = hdr;
    
end

% check if anything unequal --> error
if numel(unique(sample_count_data)) > 1
    error('Data sizes unequal.');
end

if numel(unique(sample_count_tvec)) > 1
    error('tvec sizes unequal.');
end


% % check if ExpKeys available (EC removed)
% keys_f = FindFiles('*keys.m');
% if ~isempty(keys_f)
%     run(keys_f{1});
%     csc_tsd.cfg.ExpKeys = ExpKeys;
% end

% add sessionID
[~,csc_tsd.cfg.SessionID,~] = fileparts(pwd);

% put the decimation factor back in if it was reactive instead of hard
% coded in the decimateByFactor.
if ~isempty(cfg.desired_sampling_frequency) && isempty(cfg.decimateByFactor)
    cfg.decimateByFactor = this_decimate;
end

% housekeeping
csc_tsd.cfg.history.mfun = cat(1,csc_tsd.cfg.history.mfun,mfun);
csc_tsd.cfg.history.cfg = cat(1,csc_tsd.cfg.history.cfg,{cfg});


%
function csc_info = readCSCHeader(Header)

csc_info = [];
for hline = 1:length(Header)
    
    line = strtrim(Header{hline});
    
    if isempty(line) || ~strcmp(line(1),'-') % not an informative line, skip
        continue;
    end
    
    % if we are reading an older version header with colons, remove them
    colon_idx = strfind(line,':');
    if ~isempty(colon_idx)
        line(colon_idx) = [];
    end
    
    % This expression captures the first chunk of alphanumeric characters and puts it into
    % <key> then puts whatever is to the right of it into <val>. If there is only one
    % character chunk (e.g., missing value) then it returns <val> as empty.
    a = regexp(line(2:end),'(?<key>^\S+)\s+(?<val>.*)|(?<key>\S+)','names');
    
    % deal with characters not allowed by MATLAB struct
    if strcmp(a.key,'DspFilterDelay_�s') || strcmp(a.key,'DspFilterDelay_�s') || strcmp(a.key,'DspFilterDelay_�s')
        a.key = 'DspFilterDelay_us';
    end
    
    csc_info = setfield(csc_info,a.key,a.val);
    
    % convert to double if possible
    if ~isnan(str2double(a.val))
        csc_info = setfield(csc_info,a.key,str2double(a.val));
    end
    
end