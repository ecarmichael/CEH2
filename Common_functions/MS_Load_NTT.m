function [S] = MS_Load_NTT(cfg_in, fname)
%% MS_Load_NTT: loads and extracts spike times for all cells from a .ntt file
%
%
%
%    Inputs: 
%    - cfg [struct]   configuration see the defaults below. 
%
%    - 
%
%
%
%    Outputs: 
%    -
%
%
%
%
% EC 2025-06-10   initial version 
%
%
%
%% initialize
cfg_def = [];



cfg = ProcessConfig(cfg_def, cfg_in);

if nargin < 2
    fname = dir('TT*.ntt');
    
    % remove files with 
    k_idx = false(1,length(fname)); 

    for ii = 1:length(fname)
        if contains(fname(ii).name, '000')
            k_idx(ii) = true;
        end
    end
    fname(~k_idx) = []; 
elseif nargin == 2  && (sum(contains(fname, '*')) > 0) % deal with wildcard calls.
    
    fname = dir(fname); 

end

%% loop over files and 
S = ts(); 
% S.type = 'ts';
% S.t = {}; 
% S.label = {};  
S.wave = {}; 
% S.cfg.history = {}; 

for iF = 1:length(fname)


    [Timestamps, CellNumbers, Samples, ~] =...
        Nlx2MatSpike(fname(iF).name, [1 0 1 0 1], 1, 1, [] );

    tt_id = fname(iF).name(1:strfind(fname(iF).name, '.ntt')-1);

    c_id = unique(CellNumbers);

    % loop over cells; 
    for iC = 1:length(c_id)
        S.t{end+1} = Timestamps(CellNumbers == c_id(iC))'*10^-6;
        S.label{end+1} = [tt_id '_' num2str(c_id(iC)+1)];
        S.wave{end+1} = mean(Samples(:,:, (CellNumbers == c_id(iC))),3); 
    end


end % file loop

fprintf('<strong>%s</strong>: <strong>%.0f</strong> cells loaded\n', mfilename, length(S.t))

