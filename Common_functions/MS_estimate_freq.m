function IPI_vec = MS_estimate_freq(tvec, data)
%% MS_estimate_freq: estimates the feequency using interpeak intervals.
%
%
%
%    Inputs:
%    - tvec:  [1 x nSamples]   time vector
%
%    - data: 1 x nSamples]   data vector
%
%
%    Outputs:
%    - freq_out: [1 x nSample array]  estmated frequency per cycle based on
%    interpeak interval.
%
%
%
%
% EC 2021-06-29   initial version
%% get the phase
 Phi = angle(hilbert(data)); 


phi_peaks = [];
for ii = 1:length(Phi)-1
    if (Phi(ii)>0) && (Phi(ii+1)<=0)
        phi_peaks = [phi_peaks ii+1];
    end
end

% get the Inter Peak Interval
IPI = diff(tvec(phi_peaks));

IPI_vec = nan(size(Phi));

for ii = 1:length(phi_peaks)
    if ii == 1 % fill in the first data point to the first peak
        IPI_vec(1:phi_peaks(ii+1)) = IPI(ii);
    elseif ii == length(phi_peaks)
        IPI_vec(phi_peaks(ii):end) = IPI(ii-1);
    else
        IPI_vec(phi_peaks(ii):phi_peaks(ii+1)) = IPI(ii);
    end
end

IPI_vec = 1./IPI_vec; % convert to Hz. 

if length(IPI_vec) ~= length(tvec)
    error('IPI and original tvec differ in length');
end
