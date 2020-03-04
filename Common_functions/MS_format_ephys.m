function [csc, evt] = MS_format_ephys(cfg_in, fname)
%% MS_Load_ephys: wrapper function for converting 'ephys.mat' files to csc format for MS analyses. 
%
%
%
%    Inputs: 
%     - cfg_in [struct] contains user configuration parameters (see below)
%
%     - fname: [string] file to be converted to csc format.  Should be
%     'ephys.mat'. if empty, will look for raw .ncs and .nev files and load
%     them using the 
%
%     - 
%
%
%
%    Outputs: 
%     - csc [struct] data structure with NLX continuously sampled channels
%     ('csc') 
%
%     - evt [struct] data structure with NLX event times and labels. 
%
%
% EC 2020-02-27   initial version 
%
%
%% initialize. 

if nargin == 0
    
    
elseif nargin ==1
    waring('No cfg_in specified, using defaults...')
end


cfg_def = [];
cfg_def.fc = {};
cfg_def.label = {};
cfg_def.TimeConvFactor = 10^-6; % 10^-6 means convert nlx units to seconds
cfg_def.VoltageConvFactor = 1; % 1 means output in volts, 1000 in mV, 10^6 in uV
cfg_def.decimateByFactor = []; % how much would you like to
cfg_def.desired_sampling_frequency = [];  % put in what you want the final FS to be.  Commonly 2kHz


cfg = ProcessConfig(cfg_def, cfg_in);


%% find and load the data;

if exist(fname, 'file')
    fprintf('<strong>MS_format_ephys:</strong> loading <strong>%s</strong>...\n', fname)
    this_data = load(fname);

else
    
    
end

%get the Nlx2MatCSC vars
Timestamps = this_data.data.LFP.Nlx2MatCSC.TimeStamps;
Samples = this_data.data.LFP.Nlx2MatCSC.Samples; 

NumberOfValidSamples = this_data.data.LFP.Nlx2MatCSC.NumberOfValidSamples; 

%% get the minimal NLX header information

    % extract information from header
    hdr = readCSCHeader(this_data.data.LFP.Nlx2MatCSC.CSHeader);
    
    Fs = hdr.SamplingFrequency;


%% convert the timestamps from blocks to a full tvec



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
    
   fprintf('<strong>%s</strong>: %s %d/%d bad blocks found (%.2f%%).\n',mfilename,fname,badBlocks,nBlocks,(badBlocks./nBlocks).*100);
    
    
    % remove nans
    data = data(~isnan(data));
    tvec = tvec(~isnan(tvec));

    
    
    %% check for downsammpling and decimate if needed. 
    
    if this_data.data.LFP.SamplingRate ~= Fs
        fprintf('<strong>%s</strong> Matching downsampling', mfilename)
        
        Fs_diff = Fs /this_data.data.LFP.SamplingRate; 
        
        tvec_down = decimate(tvec, Fs_diff);
        
        
        
        
%% convert to 'tsd' format for the CSC


for iChan = length(this_data.data.LFP.ChannelNames):-1:1
    
    

%     csc = tsd(this_data.data.

end


%% convert to 'ts' format for evt. 





% functions
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
    if strcmp(a.key,'DspFilterDelay_�s') || strcmp(a.key,'DspFilterDelay_�s') || strcmp(a.key, 'DspFilterDelay_µs')
        a.key = 'DspFilterDelay_us';
    end
    
    csc_info = setfield(csc_info,a.key,a.val);
    
    % convert to double if possible
    if ~isnan(str2double(a.val))
        csc_info = setfield(csc_info,a.key,str2double(a.val));
    end
    
end
