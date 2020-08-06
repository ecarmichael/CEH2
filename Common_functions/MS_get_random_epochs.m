function IV_out = MS_get_random_epochs(tvec, varargin)
%% MS_get_random_epochs: Generate a set of non-overlapping duration matched ephochs
%
%
%
%    Inputs:
%   - tvec [n x samples] time vector to draw IV's from
%
%   % if only two inputs
%   - IV [struct] IV (interval) data struct containing start and stop times.
%
%   % if 3 inputs
%   - tstart [array]  list of start times
%
%   - tend: [array] list of end times
%
%
%   Outputs:
%    - IV_out: [struct] contains IV formated start and stop times.
%
%
%
%
% EC 2020-05-08   initial version
%
%
%
%% initialize

if nargin == 2
    if ~strcmp(varargin{1}.type,'iv')
        error('IV input is not in the ''IV'' format. use iv function');
    end
    
    IV_in= varargin{1};
    
elseif nargin == 3
    
    IV_in = iv(varargin{1},varargin{2});
    
else
    error('Requires tvec and either IV or tstart & tend variables. ');
end

%% generate some random intervals.

% get idx for tstart and tend
S_idx = ismember(tvec, IV_in.tstart');
E_idx =  ismember(tvec, IV_in.tend');

success = 0;
while success == 0 % keep shifting until non of the start and stops cross the shifted end and
    %randomly shift tvec
    shift_tvec = circshift(tvec, floor(MS_randn_range(1,1,1,length(tvec))));
    
    
    % rebuild the IV using S_idx and E_idx to preserve duration.
    % sort the values. once you get it right.
    [tend, idx] = sort(shift_tvec(E_idx));
    tstart = shift_tvec(S_idx);
    tstart = tstart(idx); 

    if min(tend - tstart) >=0 % check that nothing is strange around the shifted end of the tvec.
        IV_out = iv(tstart, tend);
        
        CheckIV(IV_out);
        success = 1;
        return
    end
    
end
