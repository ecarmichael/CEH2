function ModIdx = MS_ModIdx(phi_amp)
%% MS_ModIdx: compute the modulation index for a Phase - amplitude distribution 'phi_amp'
%
%
%
%    Inputs: 
%    - phi_amp: [1 x nSamples] the mean amplitude x phase bins
%
%
%
%    Outputs: 
%    - ModIdx: [double] modulation index from Tort et al. 
%
%
%
%
% EC 2022-04-27   initial version 
%
%
%
%% compute Mod index

nbins = length(phi_amp); 

ModIdx = (log(nbins)-(-sum((phi_amp/sum(phi_amp)).*log((phi_amp/sum(phi_amp))))))/log(nbins); 
