function SWR_align_Ca_NLX(ca_dir, nlx_dir, save_dir)
%% SWR_align_Ca_NLX:
%
%
%
%    Inputs: 
%    - ca_dir: [path]      path to dir containing the ms files to align. if
%    empty will look in current directory
%
%    - nlx_dir: [path]     path to nlx data. if empty will look in current
%    dir
%
%    - save_dir: [path]  path to save the output. 
%
%    Outputs: 
%    - none
%
%
%
%
% EC 2023-10-21   initial version 
%
%
%
%% Get the nlx data

cd(nlx_dir)

evts = LoadEvents([]);

rm_idx = contains(evts.label, 'St');
evts.t(rm_idx) = [];
evts.label(rm_idx) = [];

for ii = 1:length(evts.label)
    fprintf('Event %d (%d samples)\n', ii, length(evts.t{ii}))

    
end


%% get the ca_dir time stamps
Ca = []; warning off;

for ii = length(ca_dir):-1:1
    cd(ca_dir{ii})
    Ca{ii} =   load('ms.mat', 'ms');
    parts = strsplit(ca_dir{ii}, filesep);
    parts = strsplit(parts{end}, '_');
    Ca{ii}.ms.info.subject = parts{1};
    Ca{ii}.ms.info.date = [parts{2} '_' parts{3} '_' parts{4}];
    
    parts = strsplit(ca_dir{ii}, filesep);
    pdir = parts{end-1};
    parts = strsplit(pdir, '_');
    Ca{ii}.ms.info.sess = parts{end}; 
    
    if strcmpi(Ca{ii}.ms.info.subject, '1860')
        Ca{ii}.ms.nlx.tvec = sort([unique(evts.t{1}), unique(evts.t{2})]);
        
        cd(nlx_dir)
        cfg_csc.fc = {'CSC1.ncs', 'CSC5.ncs', 'CSC8.ncs'};
        cfc_csc.desired_sampling_frequency = 2000;
        Ca{ii}.ms.nlx.csc = MS_LoadCSC(cfg_csc);
        Ca{ii}.ms.nlx.labels = {'Contra1', 'Ipsi1', 'Ipsi2'}; 
        
        
        
    elseif strcmpi(Ca{ii}.ms.info.subject, '1827')
                cd(nlx_dir)
        Ca{ii}.ms.nlx.tvec = sort([unique(evts.t{3}), unique(evts.t{4})]);
        cfg_csc.fc = {'CSC9.ncs', 'CSC11.ncs', 'CSC10.ncs'};
        cfc_csc.desired_sampling_frequency = 2000;
        Ca{ii}.ms.nlx.csc = MS_LoadCSC(cfg_csc);
        Ca{ii}.ms.nlx.labels = {'Contra1', 'Ipsi1', 'Ipsi2'}; 
        
    end
    
end

warning on

%% Check the TS from the Ca

for ii = 1:length(Ca)
    
        fprintf('<strong>%s %s</strong>  (%d samples)\n', Ca{ii}.ms.info.subject, Ca{ii}.ms.info.sess, length(Ca{ii}.ms.time))
        
end

% align the nlx to the Ca frames by filling in missing values. 

%% simplify the ms 

for ii = length(ca_dir):-1:1
    
    ms = [];
    ms.info = Ca{ii}.ms.info;
    ms.numNeurons = Ca{ii}.ms.numNeurons;
    ms.Options = Ca{ii}.ms.Options;
    ms.time = Ca{ii}.ms.time;
    ms.tvecs = Ca{ii}.ms.tvecs;
    ms.Centroids = Ca{ii}.ms.Centroids;
    ms.CorrProj = Ca{ii}.ms.CorrProj;
    ms.SFPs = Ca{ii}.ms.SFPs;
    ms.SFPs_sharp = Ca{ii}.ms.SFPs_sharp;
    ms.detrendRaw = Ca{ii}.ms.detrendRaw;
    ms.FiltTraces = Ca{ii}.ms.FiltTraces;
    ms.RawTraces = Ca{ii}.ms.RawTraces;
    ms.detrendRaw = Ca{ii}.ms.detrendRaw;
    ms.Binary = Ca{ii}.ms.Binary;
    ms.denoise = Ca{ii}.ms.denoise;
    ms.deconv = Ca{ii}.ms.deconv;
    ms.detrendRaw = Ca{ii}.ms.detrendRaw;
    ms.nlx = Ca{ii}.ms.nlx; 
   save([save_dir filesep ms.info.subject '_' ms.info.date '_' ms.info.sess '.mat'], 'ms', '-v7.3')
    
end




