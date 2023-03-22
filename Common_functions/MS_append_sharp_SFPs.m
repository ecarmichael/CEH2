function ms = MS_append_sharp_SFPs(ms)
%% MS_append_sharp_SFPs: sharpends the SFPs in the ms struct using the Ziv method. 
%
%    based on msExtractSFPs by Guillaume Etter Contact: etterguillaume@gmail.com
%
%
%
%
%    Inputs: 
%    - ms : [struct] 'ms' structure from cnmfe pipeline
%
%
%
%    Outputs: 
%    - ms: [struct] 'ms' with the 'SFPs_sharp' field added. 
%
%
%
% EC 2023-03-21   initial version 
%
%
%
%% 
ms.SFPs_sharp = [];

for ii = size(ms.SFPs,3):-1:1
    SFP_t = ms.SFPs(:,:,ii);
    SFP_t(SFP_t<0.5*max(max(SFP_t))) = 0; % This is to sharpen footprints, based on Ziv lab method
    ms.SFPs_sharp(:,:,ii) = SFP_t;
end

